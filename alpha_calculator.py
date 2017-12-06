
# coding: utf-8

# In[53]:


import pandas as pd
import numpy as np
from scipy import spatial
from scipy import interpolate
import math
import matplotlib.pyplot as plt
import matplotlib.tri as tri
plt.style.use('ggplot')


# In[99]:


# files to read data (filename) and write data (filename2) into, and mesh details (filename3)

filename = '/Users/anamolpundle/Dropbox/Fluent Data/DPM Data/200mm_particle_history_small.his'
filename2 = '/Users/anamolpundle/Dropbox/Fluent Data/DPM Data/write_here.his'
filename3 = '/Users/anamolpundle/Dropbox/Fluent Data/DPM Data/200mm_cell_id.csv'


# In[101]:


# remove header from filename and write to filename2

file = open(filename, 'r')
file2 = open(filename2, 'a')

limit = 14
i = 0
for line in file:
    if (i > limit):
        #print(line)
        file2.write(line)
    i = i + 1
    
file.close()
file2.close()


# In[115]:


# read filename2 into dataframe

df = pd.read_csv(filename2, sep = "\s+", header = None)
df.columns = ['ResTime', 'ID', 'Temperature', 'Xpos', 'Radpos']


# In[ ]:





# In[103]:


# read filename3 (mesh data) into dataframe

mesh_data = pd.read_csv(filename3)
mesh_data = mesh_data.iloc[:,0:4]
mesh_data.columns = ['cell_number', 'x-coordinate', 'y-coordinate', 'cell_id']


# In[116]:


df.head()


# In[105]:


df.shape


# In[11]:


df.describe()


# In[106]:


mesh_data.cell_id = mesh_data.cell_id.astype(int)
mesh_data.head()


# In[107]:


mesh_data.describe()


# In[109]:


def filter_func(x):
    # removes all particles that do not exit at the outlet
    goes_to_the_end = x['Xpos'].max() > 0.75
    return (goes_to_the_end)   

def calculate_thermal_age(x):
    # numerical integration of Temperature wrt residence time to calculate thermal age
    over_nucleation_zone = False
    x = x.sort_values('ResTime')
    x['thermal_age'] = 0.0
    rowsize = x.shape[0]
    for i in range(1,rowsize):
        if (x.iloc[i,3] > 0.15):
            over_nucleation_zone = True
        
        if (over_nucleation_zone):
            x.iloc[i,5] = x.iloc[i-1,5] + 0.5 * (x.iloc[i,2] + x.iloc[i-1,2])*(x.iloc[i,0] - x.iloc[i-1,0])
    return x


# In[110]:


# apply filter_func and calculate_thermal_age to dataframe after grouping by ID

dfgp = df.groupby('ID').filter(filter_func).groupby('ID').apply(calculate_thermal_age)


# In[111]:


#dfgp = dfgp.sort_values('ID')
dfgp = dfgp.reset_index(drop=True)
dfgp.head()


# In[112]:


dfgp.shape


# In[113]:


dfgp.describe()


# In[114]:


# find cell values for each point in the DPM dataset, by figuring out the nearest cell centroid
# neighbor from the mesh_data dataframe

mesh_x_and_y = np.asarray(mesh_data.iloc[:,1:3])
points_to_look_up = np.asarray(dfgp.iloc[:,3:5])

mesh_tree = spatial.KDTree(mesh_x_and_y)

cell_ids = mesh_tree.query(points_to_look_up)


# In[80]:


cell_id_array = np.asarray(cell_ids)[1,:]
cell_id_array = cell_id_array.astype(int)
a = mesh_data['cell_number'][cell_id_array]
a = a.reset_index(drop=True)
#dfgp.shape
#dfgp['cell_id'] = pd.Series(a,index=dfgp.index)
#dfgp.head()

dfgp['cell_id'] = a


# In[81]:


dfgp.head()


# In[82]:


dfgp.describe()


# In[83]:


def calc_std_and_mean(x):
    # calculates the mean and std deviation of the thermal age for a cell
    mean_thermal_age = x['thermal_age'].mean()
    std_thermal_age = x['thermal_age'].std()
    cell_id = x.iloc[0,6]
    dictionary = {'cell_id': [cell_id], 'mean_thermal_age':[mean_thermal_age], 'std_thermal_age':[std_thermal_age]}
    new_df = pd.DataFrame(data = dictionary)
    return new_df


# In[84]:


cell_values = dfgp.groupby('cell_id').apply(calc_std_and_mean)


# In[85]:


cell_values['std_thermal_age'] = cell_values['std_thermal_age'].fillna(0.0)


# In[86]:


cell_values.columns = ['cell_number', 'mean_thermal_age', 'std_thermal_age']
cell_values.describe()


# In[ ]:





# In[87]:


cell_values_ultimate = pd.merge(mesh_data, cell_values, on = 'cell_number', how='left')
cell_values_ultimate = cell_values_ultimate.drop('cell_id', axis = 1)
cell_values_ultimate = cell_values_ultimate.drop('x-coordinate', axis = 1)
cell_values_ultimate = cell_values_ultimate.drop('y-coordinate', axis = 1)


# In[88]:


cell_values_ultimate.describe()


# In[89]:


#find cell id of point with maximum soot

abc = mesh_data.loc[mesh_data['x-coordinate'] < 0.22]
pqr = abc.loc[abc['x-coordinate'] > 0.21]
xyz = pqr.loc[pqr['y-coordinate'] < 0.0004]
xyz


# In[90]:


cell_values_ultimate.loc[cell_values_ultimate['cell_number'] == 37886]


# In[91]:


# calculate alpha according to formula in 'Detailed modelling of soot oxidation by O2 and OH in 
# laminar diffusion flames'; Khosousi,A. Dworkin, S.B., Proceedings of the Combustion Institute 

max_thermal_age = 11.37
cell_values_ultimate['alpha'] = max_thermal_age * max_thermal_age /(cell_values_ultimate['mean_thermal_age'] ** 2) * (np.e ** (2 * (1 - max_thermal_age/cell_values_ultimate['mean_thermal_age'])))


# In[92]:


# fill resulting NaN alpha, mean and std values with 0.001

cell_values_ultimate['alpha'] = cell_values_ultimate['alpha'].fillna(0.001)
cell_values_ultimate['mean_thermal_age'] = cell_values_ultimate['mean_thermal_age'].fillna(0.001)
cell_values_ultimate['std_thermal_age'] = cell_values_ultimate['std_thermal_age'].fillna(0.001)


# In[93]:


cell_values_ultimate.head()


# In[94]:


cell_values_ultimate.describe()


# In[95]:


newdf = pd.merge(cell_values_ultimate, mesh_data, on='cell_number', how='outer')


# In[96]:


newdf = newdf.round(6)
newdf.describe()


# In[97]:


newdf = newdf.drop('cell_id', axis =1)
newdf.head()


# In[61]:


newdf.to_csv('/Users/anamolpundle/Dropbox/Fluent Data/DPM Data/300mm_cell_values_alpha.csv', header=False, index=False)


# In[ ]:





# In[ ]:





# In[98]:


#f, ax = plt.subplots(1,2, sharex=True, sharey=True)
error_percentage = newdf['std_thermal_age']/newdf['mean_thermal_age']
plt.style.use('ggplot')
x = newdf['y-coordinate']
y = newdf['x-coordinate']
z = newdf['alpha']
#z = error_percentage
x1 = -x
y1 = y
x_new = np.hstack([x,x1])
y_new = np.hstack([y,y1])
z_new = np.hstack([z,z])



"""ax[0].tripcolor(x,y,z)
ax[0].set_xlim[0,0.2]
#plt.colorbar()
ax[1].tricontourf(x,y,z, 100) # choose 20 contour levels, just to show how good its interpolation is
#ax[1].plot(x,y, 'ko ')
#ax[0].plot(x,y, 'ko ')
plt.show()"""
fig = plt.figure(figsize=(12,16))
plt.tripcolor(x_new,y_new,z_new)
#plt.ylabel('ylabel', fontsize=16)
#plt.tricontourf(x,y,z, 100)
plt.colorbar()
plt.show()


# In[ ]:




