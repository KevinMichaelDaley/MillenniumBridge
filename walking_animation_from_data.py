import math, mathutils
import numpy as np
import bpy
def dist(v0, v1): 
    return (mathutils.Vector(v0)-mathutils.Vector(v1)).length

class Dummy:
    def __init__(self, t, t_next_step, displace_com, U,  armature):
        self.arm=armature
        self.numstring=''
        if '.' in self.arm.name:
            self.numstring='.'+armature.name.split('.')[1]
        self.ob=bpy.data.objects['Dummy'+self.numstring]
        self.pose=armature.pose
        self.displace_com=displace_com
        self.left_leg_is_stance_leg=self.displace_com[0]<0
        self.time_index = 0
        self.t = t
        self.t_next_step = t_next_step
        self.displace_foot = U
        self.t0 = 0
        self.t1 = np.where(np.diff(self.t_next_step))[0][0]
        self.U=bpy.data.objects['U'+self.numstring]
        self.V=bpy.data.objects['V'+self.numstring]
        self.Y=bpy.data.objects['Y'+self.numstring]

        self.min_foot_offset=np.min(np.abs(U-displace_com))
        self.u0 = self.min_foot_offset
        self.u00 = self.min_foot_offset
        self.u1=self.displace_foot[0]
        self.l0= self.U.location[0:3]
        self.r0= self.V.location[0:3]
        
        dmin=min(np.abs(self.U.location[0]-self.Y.location[0]),np.abs(self.V.location[0]-self.Y.location[0]))
        self.Y.location[0]=(self.U.location[0]+self.V.location[0])/2.0
        self.U.location[0]=-1.02*dmin+self.Y.location[0]
        self.V.location[0]=1.02*(dmin)+self.Y.location[0]
        self.V.location[1]=self.U.location[1]=self.Y.location[1]
        
        self.l0= self.U.location[0:3]
        self.r0= self.V.location[0:3]
        self.s0= self.Y.location[0:3]
        self.steps=0    
        self.scale=np.abs(self.s0[2]-self.l0[2])  
    def step(self):
        if np.sign(self.displace_com[self.time_index])!=np.sign(self.displace_com[self.time_index-1]): #if we've hit the next switching time
            
            self.t0=self.time_index 
            try:
                self.t1=np.where(np.diff(np.signbit(self.displace_com[self.time_index:])))[0][0]
            except:
                return 0
            
            if self.left_leg_is_stance_leg:
                self.steps+=1
                self.V.keyframe_insert(data_path="location", frame=self.t0)
                self.V.keyframe_insert(data_path="location", frame=self.t1)
                self.U.keyframe_insert(data_path="location", frame=self.t0)
                self.l0=[self.l0[0],self.l0[1]-0.05*self.scale,self.l0[2]]
                self.U.location=[self.l0[0]-0.12*self.scale, self.l0[1], self.l0[2]]
                self.U.keyframe_insert(data_path="location", frame=self.t1)
            else:
                self.steps+=1
                self.U.keyframe_insert(data_path="location", frame=self.t0)
                self.U.keyframe_insert(data_path="location", frame=self.t1)
                self.V.keyframe_insert(data_path="location", frame=self.t0)
                self.r0=[self.r0[0],self.r0[1]-0.05*self.scale,self.r0[2]]
                self.V.location=[self.r0[0]+0.12*self.scale, self.r0[1], self.r0[2]]
                self.V.keyframe_insert(data_path="location", frame=self.t1)
            self.ob.keyframe_insert(data_path="location", frame=self.t0)
            self.arm.keyframe_insert(data_path="location", frame=self.t0)
            self.ob.location[1]-=0.05/2.0*self.scale
            self.arm.location[1]-=0.05/2.0*self.scale
            self.s0=[self.s0[0],self.s0[1]-0.05/2.0*self.scale,self.s0[2]]
            self.ob.keyframe_insert(data_path="location", frame=self.t1)
            self.arm.keyframe_insert(data_path="location", frame=self.t1)
            self.left_leg_is_stance_leg = False if self.left_leg_is_stance_leg else True #change the free leg to the stance leg and vice-versa
        
        self.Y.location=[self.s0[0]+ 0.15*(self.displace_com[self.time_index])*self.scale, self.s0[1], self.s0[2]] 
        self.Y.keyframe_insert(   data_path="location",  frame=self.time_index)
        self.time_index+=1
        return 1


if __name__=="__main__":
    data=np.loadtxt('/home/kmd/Dropbox/research/Dynamics/bridge/video_assets/txtdata/our_model_trace.mcol')[::8,:]   
    dummy_start_times=[0]*len(bpy.context.selected_objects)
    for arm in bpy.context.selected_objects:
        i=2
        if '.' in arm.name:
            i=int(arm.name.split('.')[1],10)   
        dummy=None
        for st in range(1, data.shape[0]):
            if i+1>data[st,0]:
                
                dummy_start_times[i]=st
            else:
                if dummy is None:
                 dummy=Dummy(np.arange(data.shape[0] )/60.0, np.diff(np.signbit(data[:,i+2]).astype(np.float32)), data[:, i+2], np.ones(data.shape[0])*2.0,  arm)
                if dummy.step() == 0:
                 break
        obs=[dummy.U, dummy.V, dummy.Y, dummy.ob, dummy.arm]

        T=dummy_start_times[i]
        h_fall=100
        t_fall=math.sqrt(h_fall/4.9)
        zloc=[]
        for ob in obs:
            zloc.append(ob.location[2])
        interval=data[1,1]-data[0,1]
        for st in range(int(T-t_fall/interval),T):
            h=max(0,h_fall-max(0,4.9*(st-t_fall)**2))
            for j,ob in enumerate(obs):
                ob.location[2]=zloc[j]+h
                ob.keyframe_insert(data_path="location", frame=st, index=2)
                       
    