import initsgl 
import initsfwd 
import initsadj 
from numpy import *
from scipy import interpolate
class lpTransform:
	def __init__(self,N,Nproj,Nslices,filter_type,pad):
		self.N=N
		self.Nslices=Nslices
		self.filter_type=filter_type	
		#size after zero padding in the angle direction (for nondense sampling rate)
		self.osangles=int(max(round(3.0*N/2.0/Nproj),1))
		self.Nproj=self.osangles*Nproj
		#size after zero padding in radial direction
		if (pad): 
			self.Npad=3*N/2
		else: 
			self.Npad=N

	def precompute(self):		
		#precompute parameters for the lp method
		self.Pgl=initsgl.create_gl(self.Npad,self.Nproj)
		self.Pfwd=initsfwd.create_fwd(self.Pgl)
		self.Padj=initsadj.create_adj(self.Pgl,self.filter_type)

	def fwd(self,f):
		N=self.Pgl.N
		Nproj=self.Pgl.Nproj
		Ntheta=self.Pgl.Ntheta
		Nrho=self.Pgl.Nrho
		Nslices=f.shape[0]
		N1=f.shape[1] # N1==N for no padding

		df=zeros([Nslices,N,N],dtype='float32')
		dR=zeros([Nslices,N,Nproj],dtype='float32')

		df[:,N/2-N1/2:N/2+N1/2,N/2-N1/2:N/2+N1/2]=f				
		rho=transpose(tile(self.Pgl.rhosp,[Ntheta,1]))
		for islice in range(0,Nslices): 
			flp0=zeros([Nrho,Ntheta],dtype='float32')
			R0=zeros([N,Nproj],dtype='float32')
			interp_c=interpolate.RectBivariateSpline(arange(0,N), arange(0,N), df[islice,:,:], kx=3,ky=3)		
			for k in range(0,3):	
				flp0[unravel_index(self.Pfwd.cids,(Nrho,Ntheta))]=interp_c(self.Pfwd.lp2C2[k],self.Pfwd.lp2C1[k],grid=False)
				Rlp=fft.irfft2(fft.rfft2(flp0*exp(rho))*self.Pfwd.fZ)
				interp_lp=interpolate.RectBivariateSpline(arange(0,Nrho), arange(0,Ntheta), Rlp, kx=3,ky=3)		
				R0[unravel_index(self.Pfwd.pids[k],(N,Nproj))]=interp_lp(self.Pfwd.p2lp2[k],self.Pfwd.p2lp1[k],grid=False)	
			dR[islice,:,:]=R0
		R=dR[:,N/2-N1/2:N/2+N1/2,0::self.osangles]
		return R

	def adj(self,R,cor):
		N=self.Pgl.N
		Nproj=self.Pgl.Nproj
		Ntheta=self.Pgl.Ntheta
		Nrho=self.Pgl.Nrho
		Nslices=R.shape[0]
		N1=R.shape[1] # N1==N for no padding

		df=zeros([Nslices,N,N],dtype='float32')
		dR=zeros([Nslices,N,Nproj],dtype='float32')

		shift=cor-N1/2
		dR[:,N/2-N1/2+shift:N/2+N1/2+shift,0::self.osangles]=R		
		#padding
		dR[:,:N/2-N1/2+shift,:]=tile(reshape(dR[:,N/2-N1/2+shift,:],[Nslices,1,Nproj]),[1,N/2-N1/2+shift,1])
		dR[:,N/2+N1/2+shift:,:]=tile(reshape(dR[:,N/2+N1/2+shift-1,:],[Nslices,1,Nproj]),[1,N/2-N1/2-shift,1])
	
		if (self.Padj.filter is not None):
			dRos=zeros([Nslices,4*N,Nproj],dtype='float32')
			dRos[:,2*N-N/2:2*N+N/2,:]=dR
			dRos=real(fft.fftshift(fft.ifft(fft.fftshift(fft.fftshift(fft.fft(fft.fftshift(dRos,1),axis=1),1)*tile(transpose(tile(self.Padj.filter,[Nproj,1])),[Nslices,1,1]),1),axis=1),1))
			dR=dRos[:,2*N-N/2:2*N+N/2,:]

		for islice in range(0,Nslices): 
			Rlp0=zeros([Nrho,Ntheta],dtype='float32')
			f0=zeros([N,N],dtype='float32')
			interp_p=interpolate.RectBivariateSpline(arange(0,N), arange(0,Nproj), dR[islice,:,:], kx=3,ky=3)		
			for k in range(0,3):	
				Rlp0[unravel_index(self.Padj.lpids,(Nrho,Ntheta))]=interp_p(self.Padj.lp2p2[k],self.Padj.lp2p1[k],grid=False)
				Rlp0[unravel_index(self.Padj.wids,(Nrho,Ntheta))]=interp_p(self.Padj.lp2p2w[k],self.Padj.lp2p1w[k],grid=False)
				flp=fft.irfft2(fft.rfft2(Rlp0)*self.Padj.fZ)	
				interp_lp=interpolate.RectBivariateSpline(arange(0,Nrho), arange(0,Ntheta), flp, kx=3,ky=3)		
				f0[unravel_index(self.Padj.cids,(N,N))]+=interp_lp(self.Padj.C2lp2[k],self.Padj.C2lp1[k],grid=False)	
			df[islice,:,:]=f0
		f=df[:,N/2-N1/2:N/2+N1/2,N/2-N1/2:N/2+N1/2]*self.osangles*N1/N
                return f

