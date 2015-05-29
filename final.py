"""
A simple example of an animated plot... In 3D!

see sources in Monte function
TODO
-why is there a tendency to drift right?
-is there a way to make a sort of heatmap of flux on the surface plot?
-Make parameters enterable - Tkinter?
-writeup - include flux vs. distance, flux vs. number of particles, flux vs. isotropy factor
-write about the assumptions - uses standard monte carlo but also two-term Henyey-Greenstein phase functions
-how it differs from the 741
-real attenuation absorption scattering coeff.
"""
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from matplotlib import cm
from random import randint,random

from Tkinter import *
fields = ('mu_s (cm-1)', 'mu_a (cm-1)', 'Isotropy', 'Number Particles', 'World Size (cm)','Flux Width x (cm)','Flux Width y (cm)','Flux Distance z (cm)',
          'Russian Roulette Factor','Initial E-Gamma (MeV)','Cutoff E-Gamma (MeV)','Flux v Distance Range (cm)','ALPHA')
def get_mu_a(e_gamma_0,mu_s):
    MATCHER=[
        (	1.00E-003	,	3.71E+003	),
        (	1.04E-003	,	3.39E+003	),
        (	1.07E-003	,	3.09E+003	),
        (	1.50E-003	,	1.25E+003	),
        (	2.00E-003	,	5.60E+002	),
        (	2.15E-003	,	4.58E+002	),
        (	2.30E-003	,	3.80E+002	),
        (	2.47E-003	,	3.10E+002	),
        (	2.64E-003	,	2.61E+002	),
        (	2.82E-003	,	2.16E+002	),
        (	3.00E-003	,	1.84E+002	),
        (	3.61E-003	,	1.07E+002	),
        (	4.00E-003	,	8.16E+001	),
        (	5.00E-003	,	4.22E+001	),
        (	6.00E-003	,	2.46E+001	),
        (	8.00E-003	,	1.04E+001	),
        (	1.00E-002	,	5.38E+000	),
        (	1.50E-002	,	1.70E+000	),
        (	2.00E-002	,	8.23E-001	),
        (	3.00E-002	,	3.79E-001	),
        (	4.00E-002	,	2.69E-001	),
        (	5.00E-002	,	2.26E-001	),
        (	6.00E-002	,	2.05E-001	),
        (	8.00E-002	,	1.82E-001	),
        (	1.00E-001	,	1.69E-001	),
        (	1.50E-001	,	1.49E-001	),
        (	2.00E-001	,	1.36E-001	),
        (	3.00E-001	,	1.18E-001	),
        (	4.00E-001	,	1.05E-001	),
        (	5.00E-001	,	9.60E-002	),
        (	6.00E-001	,	8.87E-002	),
        (	8.00E-001	,	7.79E-002	),
        (	1.00E+000	,	7.01E-002	),
        (	1.25E+000	,	6.27E-002	),
        (	1.50E+000	,	5.70E-002	),
        (	2.00E+000	,	4.90E-002	),
        (	3.00E+000	,	3.93E-002	),
        (	4.00E+000	,	3.37E-002	),
        (	5.00E+000	,	3.00E-002	),
        (	6.00E+000	,	2.74E-002	),
        (	8.00E+000	,	2.40E-002	),
        (	1.00E+001	,	2.19E-002	),
        (	1.50E+001	,	1.92E-002	),
        (	2.00E+001	,	1.79E-002	),
    ]
    for e,mu_s in MATCHER:
        if e==e_gamma_0:
            mu=mu_s
    return mu

def runit(entries):
   # period rate:
   mu_s = float(entries['mu_s (cm-1)'].get())
   mu_a =  float(entries['mu_a (cm-1)'].get())
   e_gamma_0 = float(entries['Initial E-Gamma (MeV)'].get())
   mu_s=get_mu_a(e_gamma_0,mu_s)
   num_particles = int(entries['Number Particles'].get())
   fluxsurface_len_x =  int(entries['Flux Width x (cm)'].get())
   fluxsurface_len_y = int(entries['Flux Width y (cm)'].get())
   flux_dist_from_src_z =  float(entries['Flux Distance z (cm)'].get())
   g = float(entries['Isotropy'].get())
   lim =  float(entries['World Size (cm)'].get())
   m = float(entries['Russian Roulette Factor'].get())
   lines_or_points=int(v0.get())
   rand_init_dir=int(v1.get())
   plot_flux_surface=int(v2.get())
   sampling=v5.get()
   cutoff_e = float(entries['Cutoff E-Gamma (MeV)'].get())
   #alpha = float(entries['ALPHA'].get())
   #print rand_init_dir,lines_or_points#mu_s,mu_a,g,num_particles,fluxsurface_len_x,fluxsurface_len_y,flux_dist_from_src_z,lim,lines_or_points

    
   fig = plt.figure()
   ax = p3.Axes3D(fig)
   
   
   ax.set_xlim3d([-lim, lim])
   ax.set_xlabel('X')
   
   ax.set_ylim3d([-lim, lim])
   ax.set_ylabel('Y')
   
   ax.set_zlim3d([-lim, lim])
   ax.set_zlabel('Z')
   
   ax.set_title('3D Monte Carlo Walking')
   
   #================
   #Some I.C.s - parameters to modify

   flux_surface_area = float(fluxsurface_len_x)*float(fluxsurface_len_y)
   #absorption and scattering coefficients

   #Following function does the bulk of the code - makes lines for the walking animation, OR determines particle flight distance and flux
   
   def Monte(length, dims=6,g=0.0,mu_a=0.03,mu_s=0.3,e_gamma_0=5.0,rand_init_dir=False, lines_or_points=True,sampling='H'):
      '''
      from http://en.wikipedia.org/wiki/Monte_Carlo_method_for_photon_transport
      AND
      'Successive order, multiple scattering of two-term Henyey-Greenstein phase functions'
      AND
      Wang, Jaques MCML Monte Carlo method for light transport 1995
      '''
      e_gamma=e_gamma_0#MeV
      mu_t = mu_a+mu_s
      lineData = np.empty((dims, length))
      #starting at origin
      lineData[:, 0] = 0.0,0.0,0.0, 0,0,0 #x,y,z,#crossings,crossing_x,crossing_y
      #INITIAL DIR
      if rand_init_dir:
         x_dir,y_dir,z_dir = randint(0,10),randint(0,10),randint(0,10)
         mag = np.sqrt(x_dir**2+y_dir**2+z_dir**2)
         #normalize
         mu_x,mu_y,mu_z = x_dir/mag,y_dir/mag,z_dir/mag
      else:
         mu_x,mu_y,mu_z = 0.0,0.0,1.0
      #START WALKING
      single_p_crossing_factor=0
      x_cross,y_cross=0,0
      for index in range(1, length):
         alpha = e_gamma/0.511
         zeta1,zeta2,zeta3,zeta4=random(),random(),random(),random()
         step_scalar = -np.log(zeta1)/mu_t
         if (sampling == 'H') or (sampling == 'B'):
             ct = (1/(2*g))*(1+g**2-((1-g**2)/(1-g+2*g*zeta2))**2) if g != 0 else 2*zeta2-1
         elif sampling =='K':
             while 1:#ang_diff_scat>zeta2: # if zeta2<ct, retain ct
                 #print ang_diff_scat
                 ct = -1 + 2*random()                                        
                 trm1 = (1 + ct**2)/2.0
                 trm2 = 1.0/(1 + alpha*(1 - ct))
                 trm3 = 2*alpha**2*(1 - ct)**2*trm2/trm1
                 ang_diff_scat = trm1*trm2**2*(1 + trm3)
                 if ang_diff_scat<random():
                     break
         phi = 2*np.pi*zeta3
         theta = np.arccos(ct)
         #print ct,theta
         c = np.sqrt(1-mu_z**2)
         st = np.sin(theta)
         sp = np.sin(phi)
         cp = np.cos(phi)
         if mu_z >= 0.999:
            mu_x,mu_y,mu_z = st*cp,st*sp,ct
         elif mu_z <= -0.999:
            mu_x,mu_y,mu_z = st*cp,-st*sp,-ct
         else:
            new_mu_x = (st/c)*(mu_x*mu_z*cp-mu_y*sp)+mu_x*ct
            new_mu_y = (st/c)*(mu_y*mu_z*cp-mu_x*sp)+mu_y*ct
            new_mu_z = -c*st*cp+mu_z*ct
            mu_x,mu_y,mu_z = new_mu_x,new_mu_y,new_mu_z
            #print mu_x,mu_y,mu_z,'----',mu_x**2+mu_y**2+mu_z**2  #(x and y 'directional cosines are >1, wtheck)
         x,y,z = mu_x*step_scalar,mu_y*step_scalar,mu_z*step_scalar
         old_z = lineData[2,index-1]
         avg_x,avg_y = lineData[0,index-1]+x/2.0,lineData[1,index-1]+y/2.0
         new_z = old_z+z
         if ((-fluxsurface_len_x/2<avg_x<fluxsurface_len_x/2) and (-fluxsurface_len_y/2<avg_y<fluxsurface_len_y/2) and ((old_z < flux_dist_from_src_z and new_z > flux_dist_from_src_z) or (old_z > flux_dist_from_src_z and new_z < flux_dist_from_src_z))):
            single_p_crossing_factor+=1*abs(1/np.dot([mu_x,mu_y,mu_z],[0,0,1])) #NOT THE ABSOLUTE VALUE, 'COMING BACK' WILL DECREASE THIS VALUE
            x_cross,y_cross=lineData[0, index-1] + x,lineData[1, index-1] + y
            #print old_z, new_z,'---',avg_x,avg_y,'---',mu_x,mu_y,mu_z#lineData[:, index]
         step = x,y,z
         #energy lost absorption
         de = (mu_a/mu_t)*e_gamma
         e_gamma=e_gamma-de
         #PHOTON TERMINATION WITH RUSSIAN ROULETTE METHOD FOR ENDING PHOTON
         #if zeta <(1/m), w=0. Else, w = initial weight
         #m=40.0
         if e_gamma<cutoff_e:
             if zeta4 < (1.0/m):
                e_gamma=m*e_gamma
                #final = lineData[:, :index]  #truncate the array
                #print 'PHOTON DEAD', lineData[:, index]
                lineData[:3, index] = lineData[:3, index-1] + step
                lineData[3, index]=single_p_crossing_factor
                lineData[4, index]=x_cross
                lineData[5, index]=y_cross
             else:
                 e_gamma=0#m*e_gamma
                 lineData[:3, index] = lineData[:3, index-1] + step
                 lineData[3, index]=single_p_crossing_factor
                 lineData[4, index]=x_cross
                 lineData[5, index]=y_cross
                 final = lineData[:, :index]
                 break
                 #final = lineData[:, :index]
                 #continue
         else:
            #e_gamma=1
            lineData[:3, index] = lineData[:3, index-1] + step
            lineData[3, index]=single_p_crossing_factor
            lineData[4, index]=x_cross
            lineData[5, index]=y_cross
            final = lineData[:, :index]
            #print lineData[:, index]
            continue
      return final if lines_or_points else final[:, -1]
      
   
   
   #ADDING THIS MYSELF - plots stopping point (absorption points) for gamma and fluz surface===================
   
   #total_crossing_factor_iso=0
   #total_crossing_factor_aniso=0
   if not lines_or_points:
      total_crossing_factor=0
      gx=[]
      gy=[]
      for j in range(num_particles):
        #ptdata = [Monte(25, 6,0.33*index,mu_a,mu_s,w,False,False) for index in range(4)]#z-dir, 4 cases
        #ptdata = [Monte(25, 6,0.99*index,mu_a,mu_s,w,False,False) for index in range(2)]#z-dir, 2 cases
        #ptdata = [Monte(25, 6,0.33*index,mu_a,mu_s,w,False,False) for index in range(4)]#rand-dir, 4 cases
        ptdata = [Monte(25, 6,g,mu_a,mu_s,e_gamma_0,rand_init_dir,lines_or_points,sampling) for index in range(1)]
        for i,d in enumerate(ptdata):
          total_crossing_factor = total_crossing_factor+d[3]
          ax.scatter(d[0],d[1],d[2])#, c=color_dict[i])
          #print 'd3 is the crossing factor',d[3]
          
      flux = total_crossing_factor/(num_particles*flux_surface_area)
      print flux, 'flux'

      if plot_flux_surface:
           x = np.arange(-fluxsurface_len_x/2, fluxsurface_len_x/2, 0.1)
           y = np.arange(-fluxsurface_len_y/2, fluxsurface_len_y/2, 0.1)
           X, Y = np.meshgrid(x, y)
           zs = np.array([flux_dist_from_src_z for x,y in zip(np.ravel(X), np.ravel(Y))])
           Z = zs.reshape(X.shape)
           ax.plot_surface(X, Y, Z)



   #The line 'walking' animation==================================
   
   def update_lines(num, dataLines, lines) :
      for line, data in zip(lines, dataLines):
         # NOTE: there is no .set_data() for 3 dim data...
         line.set_data(data[0:2, :num])
         line.set_3d_properties(data[2,:num])
      return 1
   if lines_or_points:
      # data for 4 'walks' (lines) of different anisotropy - 0,0.33,0.67,0.99
      #data = [Monte(25, 6,0.33*index,mu_a,mu_s,w,False,True) for index in range(4)]#initial z-dir
      data = [Monte(25, 6,g,mu_a,mu_s,e_gamma_0,rand_init_dir,True,sampling) for index in range(num_particles)]#rand-dir
      
      # Creating line objects.
      # NOTE: Can't pass empty arrays into 3d version of plot()
      lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for i,dat in enumerate(data)]#, label='g = %s'%(i*0.33))[0] for i,dat in enumerate(data)]
      
      # Creating the Animation object 3rd param must match Gen_RandLine(length
      #interval is how fast animation occurs
      
      line_ani = animation.FuncAnimation(fig, update_lines, 25, fargs=(data, lines),interval=50, blit=False)
   
   
   
   plt.legend()
   plt.show()



def fluxit(entries):
   # period rate:
   mu_s = float(entries['mu_s (cm-1)'].get())
   mu_a =  float(entries['mu_a (cm-1)'].get())
   e_gamma_0 = float(entries['Initial E-Gamma (MeV)'].get())
   mu_s=get_mu_a(e_gamma_0,mu_s)
   num_particles = int(entries['Number Particles'].get())
   fluxsurface_len_x =  int(entries['Flux Width x (cm)'].get())
   fluxsurface_len_y = int(entries['Flux Width y (cm)'].get())
   flux_dist_from_src_z =  float(entries['Flux Distance z (cm)'].get())
   g = float(entries['Isotropy'].get())
   lim =  float(entries['World Size (cm)'].get())
   m = float(entries['Russian Roulette Factor'].get())
   lines_or_points=int(v0.get())
   rand_init_dir=int(v1.get())
   plot_flux_surface=int(v2.get())
   sampling=v5.get()
   cutoff_e = float(entries['Cutoff E-Gamma (MeV)'].get())
   #alpha = float(entries['ALPHA'].get())
   fluxrange = int(entries['Flux v Distance Range (cm)'].get())
   #print rand_init_dir,lines_or_points#mu_s,mu_a,g,num_particles,fluxsurface_len_x,fluxsurface_len_y,flux_dist_from_src_z,lim,lines_or_points


   #fig = plt.figure()

   fig, ax = plt.subplots()
   '''
   ax = p3.Axes3D(fig)
   ax.set_xlim3d([-lim, lim])
   ax.set_xlabel('X')

   ax.set_ylim3d([-lim, lim])
   ax.set_ylabel('Y')

   ax.set_zlim3d([-lim, lim])
   ax.set_zlabel('Z')

   ax.set_title('Flux v D')
   '''
   #================
   #Some I.C.s - parameters to modify

   flux_surface_area = float(fluxsurface_len_x)*float(fluxsurface_len_y)
   #absorption and scattering coefficients

   #Following function does the bulk of the code - makes lines for the walking animation, OR determines particle flight distance and flux

   def FluxMonte(length, dims=6,g=0.0,mu_a=0.03,mu_s=0.3,e_gamma_0=5.0,rand_init_dir=False, dis=1.0, sampling='H'):
      '''
      from http://en.wikipedia.org/wiki/Monte_Carlo_method_for_photon_transport
      AND
      'Successive order, multiple scattering of two-term Henyey-Greenstein phase functions'
      AND
      Wang, Jaques MCML Monte Carlo method for light transport 1995
      '''
      e_gamma=e_gamma_0#MeV
      mu_t = mu_a+mu_s
      lineData = np.empty((dims, length))
      #starting at origin
      lineData[:, 0] = 0.0,0.0,0.0, 0,0,0 #x,y,z,#crossings,crossing_x,crossing_y
      #INITIAL DIR
      if rand_init_dir:
         x_dir,y_dir,z_dir = randint(0,10),randint(0,10),randint(0,10)
         mag = np.sqrt(x_dir**2+y_dir**2+z_dir**2)
         #normalize
         mu_x,mu_y,mu_z = x_dir/mag,y_dir/mag,z_dir/mag
      else:
         mu_x,mu_y,mu_z = 0.0,0.0,1.0
      #START WALKING
      single_p_crossing_factor=0
      x_cross,y_cross=0,0
      for index in range(1, length):
         alpha = e_gamma/0.511
         zeta1,zeta2,zeta3,zeta4=random(),random(),random(),random()
         step_scalar = -np.log(zeta1)/mu_t
         if (sampling == 'H') or (sampling == 'B'):
             ct = (1/(2*g))*(1+g**2-((1-g**2)/(1-g+2*g*zeta2))**2) if g != 0 else 2*zeta2-1
         elif sampling =='K':
             while 1:#ang_diff_scat>zeta2: # if zeta2<ct, retain ct
                 #print ang_diff_scat
                 ct = -1 + 2*random()                                        
                 trm1 = (1 + ct**2)/2.0
                 trm2 = 1.0/(1 + alpha*(1 - ct))
                 trm3 = 2*alpha**2*(1 - ct)**2*trm2/trm1
                 ang_diff_scat = trm1*trm2**2*(1 + trm3)
                 if ang_diff_scat<random():
                     break       
         phi = 2*np.pi*zeta3
         theta = np.arccos(ct)
         #print ct,theta
         c = np.sqrt(1-mu_z**2)
         st = np.sin(theta)
         sp = np.sin(phi)
         cp = np.cos(phi)
         if mu_z == 1:
            mu_x,mu_y,mu_z = st*cp,st*sp,ct
         elif mu_z == -1:
            mu_x,mu_y,mu_z = st*cp,-st*sp,-ct
         else:
            new_mu_x = (st/c)*(mu_x*mu_z*cp-mu_y*sp)+mu_x*ct
            new_mu_y = (st/c)*(mu_y*mu_z*cp-mu_x*sp)+mu_y*ct
            new_mu_z = -c*st*cp+mu_z*ct
            mu_x,mu_y,mu_z = new_mu_x,new_mu_y,new_mu_z
            #print mu_x,mu_y,mu_z,'----',mu_x**2+mu_y**2+mu_z**2  #(x and y 'directional cosines are >1, wtheck)
         x,y,z = mu_x*step_scalar,mu_y*step_scalar,mu_z*step_scalar
         old_z = lineData[2,index-1]
         avg_x,avg_y = lineData[0,index-1]+x/2.0,lineData[1,index-1]+y/2.0
         new_z = old_z+z
         if ((-fluxsurface_len_x/2<avg_x<fluxsurface_len_x/2) and (-fluxsurface_len_y/2<avg_y<fluxsurface_len_y/2) and ((old_z < dis and new_z > dis) or (old_z > dis and new_z < dis))):
            single_p_crossing_factor+=1*abs(1/np.dot([mu_x,mu_y,mu_z],[0,0,1])) #NOT THE ABSOLUTE VALUE, 'COMING BACK' WILL DECREASE THIS VALUE
            x_cross,y_cross=lineData[0, index-1] + x,lineData[1, index-1] + y
            #print old_z, new_z,'---',avg_x,avg_y,'---',mu_x,mu_y,mu_z#lineData[:, index]
         step = x,y,z
         #energy lost absorption
         de = (mu_a/mu_t)*e_gamma
         e_gamma=e_gamma-de
         #PHOTON TERMINATION WITH RUSSIAN ROULETTE METHOD FOR ENDING PHOTON
         #if zeta <(1/m), w=0. Else, w = initial weight
         #m=40.0
         if e_gamma<cutoff_e:
             if zeta4 < (1.0/m):
                e_gamma=m*e_gamma
                #final = lineData[:, :index]  #truncate the array
                #print 'PHOTON DEAD', lineData[:, index]
                lineData[:3, index] = lineData[:3, index-1] + step
                lineData[3, index]=single_p_crossing_factor
                lineData[4, index]=x_cross
                lineData[5, index]=y_cross
             else:
                 e_gamma=0#m*e_gamma
                 lineData[:3, index] = lineData[:3, index-1] + step
                 lineData[3, index]=single_p_crossing_factor
                 lineData[4, index]=x_cross
                 lineData[5, index]=y_cross
                 final = lineData[:, :index]
                 break
                 #final = lineData[:, :index]
                 #continue
         else:
            #e_gamma=1
            lineData[:3, index] = lineData[:3, index-1] + step
            lineData[3, index]=single_p_crossing_factor
            lineData[4, index]=x_cross
            lineData[5, index]=y_cross
            final = lineData[:, :index]
            #print lineData[:, index]
            continue
      return final[:, -1]



   #ADDING THIS MYSELF - plots stopping point (absorption points) for gamma and fluz surface===================
   #total_crossing_factor_iso=0
   #total_crossing_factor_aniso=0
   to_plot=[]
   if sampling == 'B':
       to_plot = ['H','K']
   if sampling == 'H':
       to_plot = ['H']
   if sampling == 'K':
       to_plot = ['K']
   for p in to_plot:
       dis_a=[]
       flux_a=[]
       for dis in range(fluxrange*4):
          dis = dis/4.0
          total_crossing_factor=0
          gx=[]
          gy=[]
          for j in range(num_particles):
            #ptdata = [Monte(25, 6,0.33*index,mu_a,mu_s,w,False,False) for index in range(4)]#z-dir, 4 cases
            #ptdata = [Monte(25, 6,0.99*index,mu_a,mu_s,w,False,False) for index in range(2)]#z-dir, 2 cases
            #ptdata = [Monte(25, 6,0.33*index,mu_a,mu_s,w,False,False) for index in range(4)]#rand-dir, 4 cases
            ptdata = [FluxMonte(25, 6,g,mu_a,mu_s,e_gamma_0,rand_init_dir,dis,p) for index in range(1)]
            for i,d in enumerate(ptdata):
              total_crossing_factor = total_crossing_factor+d[3]
              #ax.scatter(d[0],d[1],d[2])#, c=color_dict[i])
              #print 'd3 is the crossing factor',d[3]

          flux = total_crossing_factor/(num_particles*flux_surface_area)
          #print 'dis',dis,'flux', flux
          dis_a.append(dis)
          flux_a.append(flux)
       plt.plot(dis_a,flux_a,label=p)
   ax.set_xlabel('Distance (cm)')
   ax.set_ylabel('Scalar Flux')
   ax.set_title('Scalar Flux v Distance')
   plt.legend()
   plt.show()








def makeform(root, fields):
   entries = {}
   for field in fields:
      #print field
      row = Frame(root)
      lab = Label(row, width=22, text=field+": ", anchor='w')
      ent = Entry(row)
      if field == 'mu_s (cm-1)':
          ent.insert(0,"0.16")
      elif field =='mu_a (cm-1)':
          ent.insert(0,"0.38")
      elif field == 'Number Particles':
          ent.insert(0,"500")
      elif field == 'Isotropy':
          ent.insert(0,"0.8")
      elif field == 'World Size (cm)':
          ent.insert(0,"40")
      elif field == 'Flux Width x (cm)':
          ent.insert(0,"1")
      elif field == 'Flux Width y (cm)':
          ent.insert(0,"1")
      elif field == 'Flux Distance z (cm)':
          ent.insert(0,"4")
      elif field == 'Russian Roulette Factor':
          ent.insert(0,"1.4")
      elif field =='Cutoff E-Gamma (MeV)':
          ent.insert(0,"0.03")
      elif field == 'Initial E-Gamma (MeV)':
          ent.insert(0,"1.25")
      elif field == 'Flux v Distance Range (cm)':
          ent.insert(0,"5")
      elif field == 'ALPHA':
          ent.insert(0,"2.0")
      row.pack(side=TOP, fill=X, padx=5, pady=5)
      lab.pack(side=LEFT)
      ent.pack(side=RIGHT, expand=YES, fill=X)
      entries[field] = ent
   return entries

if __name__ == '__main__':
   root = Tk()
   ents = makeform(root, fields)
   #walk or flux
   MODES0 = [
      ("Walking Animation", "1"),
      ("Flux Scatter Plot", "0"),
   ]
   v0 = StringVar()
   v0.set("0") # initialize
   for text, mode in MODES0:
      b0 = Radiobutton(root, text=text,
                  variable=v0, value=mode)
      b0.pack(anchor=W)
   #initial direction
   MODES1 = [
      ("Random Init Dir", "1"),
      ("Z-Directed", "0"),
   ]
   v1 = StringVar()
   v1.set("0") # initialize
   for text, mode in MODES1:
      b1 = Radiobutton(root, text=text,
                  variable=v1, value=mode)
      b1.pack(anchor=E)
   #plot flux surface??
   MODES2 = [
      ("Plot Flux Surface", "1"),
      ("Don't Plot Flux Surface", "0"),
   ]
   v2 = StringVar()
   v2.set("0") # initialize
   for text, mode in MODES2:
      b2 = Radiobutton(root, text=text,
                  variable=v2, value=mode)
      b2.pack(anchor=W)
   #plot flux surface??
   MODES5 = [
      ("Klein Nishina", "K"),
      ("Henyey Greenstein", "H"),
      ("Both (Compare - only for flux v distance)", "B"),
   ]
   v5 = StringVar()
   v5.set("H") # initialize
   for text, mode in MODES5:
      b5 = Radiobutton(root, text=text,
                  variable=v5, value=mode)
      b5.pack(anchor=E)
   #root.bind('<Flux v D>', (lambda event, e=ents: fetch(e)))
   root.bind('<Return>', (lambda event, e=ents: fetch(e)))   
   b2 = Button(root, text='3d Visual',
          command=(lambda e=ents: runit(e)))
   b2.pack(side=LEFT, padx=5, pady=5)
   b4 = Button(root, text='Flux v Distance',
          command=(lambda e=ents: fluxit(e)))
   b4.pack(side=LEFT, padx=5, pady=5)
   b3 = Button(root, text='Quit', command=root.quit)
   b3.pack(side=LEFT, padx=5, pady=5)
   #Tkinter radio buttons
   root.mainloop()


