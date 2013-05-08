from sys import *
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.widgets import *
from matplotlib.patches import *
from scipy.integrate import *
import scipy.optimize as opt


ax = plt.subplot(111)
plt.subplots_adjust(bottom=0.2)

def vmax(r,R,V0,nu):
	return V0*((r/R)**nu)

def prob_dist(f,min,max):
	fmin = 0.0
	fmax = f(0.0)
	trial1 = random.uniform(min,max,1000)
	trial2 = random.uniform(fmin,fmax,1000)
	for u,v in zip(trial1,trial2):
		if v <= f(u):
			value = u
			break
	return value
	
def Vel_field(sigcl):
	f_v = lambda v: 1.0/(sigcl*sqrt(2.*pi)) * exp(-(v**2)/(2.*sigcl**2))
	v_mag = prob_dist(f_v,-sigcl*5,sigcl*5)
	theta = random.uniform(0.0,pi)
	phi = random.uniform(0.0,2*pi) 
	return (v_mag,theta,phi)	

class main:
	def __init__(self):
		self.cmap = plt.cm.gist_rainbow
		self.R=10.
		self.Rcl = .01
		self.Rb = .05
		self.V0=1.
		self.f=.05
		self.n0 = 100 
		self.sigth = 0.1
		self.Gamma = 1.0
		self.gamma = 0.5
		self.nu=.3
		self.seednum=0
		self.N = 100
		self.field = 1
		self.sigcl = 1.0
		self.gs = plt.GridSpec(28,15)
		self.gs.update(left=0.10,right=.93, top=0.97,bottom=.21,hspace=2)
		
	def Rcl_recalc(self,Rcl):
		self.Rcl = Rcl

	def Rb_recalc(self,Rb):
		self.Rb = Rb

	def Vinit_recalc(self,Vinit):
		self.V0 = Vinit

	def n_recalc(self,n):
		self.n0 = n

	def f_recalc(self,f):
		self.f = f

	def sigcl_recalc(self,sigcl):
		self.sigcl=sigcl

	def Gamma_recalc(self,Gamma):
		self.Gamma = Gamma

	def R_recalc(self,R):
		self.R = R

	def variables(self):
		self.N = (self.R/self.Rcl)**3 * self.f
		
	def distribute(self):
		N1d = int(self.N**(1./3.))
		scale = self.N**(1./3.)/(2.*self.R)
		Nb1d = self.Rb*scale
		self.seednum += 1
		random.seed(self.seednum)		

		x_const_y,x_const_z,y_const_z,z_const_y = [],[],[],[]
		params = []
		const_y_Us = []
		const_y_Vs = []
		const_z_Us = []
		const_z_Vs = []

		field = '2'
		los = [[],[]]
		prof_vel = []
		for z in range(N1d):
			z+=1
			for y in range(-int(Nb1d)-1,int(Nb1d)+2):
				for x in range(-int(Nb1d)-1,int(Nb1d)+2):
					if y == 0 and  z <= int(Nb1d)*15+10:
						params = Vel_field(self.sigcl)
						Vx = params[0]*sin(params[1])*cos(params[2])
						Vy = params[0]*sin(params[1])*sin(params[2])
						Vz = params[0]*cos(params[1])	
						z_const_y.append(z/scale)
						x_const_y.append(x/scale)
						const_y_Us.append(Vz)
						const_y_Vs.append(Vx)
					if z == 1:
						params = Vel_field(self.sigcl)
						Vx = params[0]*sin(params[1])*cos(params[2])
						Vy = params[0]*sin(params[1])*sin(params[2])
						Vz = params[0]*cos(params[1])	
						x_const_z.append(x/scale)
						y_const_z.append(y/scale)
						const_z_Us.append(Vx)
						const_z_Vs.append(Vy)
			
		p1 = plt.subplot(self.gs[:11,2:7])
		plt.xlim([-self.Rb,self.Rb])
		plt.ylim([-self.Rb,self.Rb])
		plt.xlabel('X (pc)',fontsize=15)
		plt.ylabel('Y (pc)',fontsize=15)
		plt.quiver(x_const_z,y_const_z,const_z_Us,const_z_Vs,color = 'r', linewidths = (2,), edgecolors = ('k'), headaxislength=5, scale = 10)
		for x_cen,y_cen in zip(x_const_z,y_const_z):
			circle = Circle((x_cen,y_cen),radius=self.Rcl,alpha=0.2)
			p1.add_patch(circle)
		circle = Circle((0,0),self.Rb,fill=False,color='black')
		p1.add_patch(circle)
		p2 = plt.subplot(self.gs[13:24,:9])
		plt.xlim([0,self.Rb/2 * 10.])
		plt.ylim([-self.Rb,self.Rb])
		plt.xlabel('Z (pc)',fontsize=15)
		plt.ylabel('X (pc)',fontsize=15)
		plt.quiver([foo-.15 for foo in z_const_y],x_const_y,const_y_Us,const_y_Vs,color = 'r', linewidths = (2,), edgecolors = ('k'), headaxislength=5, scale = 20)
		for x_cen,z_cen in zip(x_const_y,z_const_y):
			circle = Circle((z_cen-.15,x_cen),self.Rcl,alpha=0.2)
			p2.add_patch(circle)
		plt.draw()
		
	def line_prof(self):
		p3 = plt.subplot(self.gs[2:22,11:])
		Nb = self.N * (pi * self.Rb**2 * 2)/((4./3.)* pi * self.R**2) 
		clump_vs = []
		for clump in range(int(Nb)):
			params = Vel_field(self.sigcl) 
			Vz = params[0]*cos(params[1])	
			clump_vs.append(Vz)
		pc = 3.08567758e18
		R_19 = self.R * pc * 1e-19
		Rcl_19 = self.Rcl * pc * 1e-19 
		Tex = 8.
		Tbb = 2.725
		T = 100.
		b = 12.9 * sqrt(T*1e-4/1.008)
		#empty_rad = (3./(4.*pi) * pi * self.Rb**2 * 2. * self.R * (1.-self.f))**(1./3.) 
		empty_rad = (3./(4.*pi) *  Nb * (4.*pi/3.) * self.Rcl**3 * ((1./self.f)-1.))**(1./3.)
		empty_rad_19 = empty_rad * pc * 1e-19		
		tau_list = []
		vbar_list = []
		Tb_list = []
		n0  = self.n0
		f = self.f
		Gamma = self.Gamma
		vbar_list = linspace(min(clump_vs)-3,max(clump_vs)+3,200)
		delt_v = ((max(clump_vs)+3) - (min(clump_vs)-3))/200.
		for v in linspace(min(clump_vs)-3,max(clump_vs)+3,200):
			tau_eff = 0
			ind_tau_list = map(lambda Vz: 3.53 * Gamma * n0 * 1e-3 * R_19 * (2.0/b) * 1./(b*sqrt(pi)) * exp(-(Vz-v)**2/(b**2))*f*(1./int(Nb)),clump_vs)
			tau_eff = sum(ind_tau_list)
			#for Vz in clump_vs:
			#	tau_eff += 3.53 * self.n0 * 1e-3 * R_19 * (2.0/b) * 1./(b*sqrt(pi)) * exp(-(Vz-v)**2/(b**2))*self.f*(1./int(Nb)) 
			tau_eff += 3.53 * n0 * 1e-3 * R_19 * (2.0/b) * 1./(b*sqrt(pi)) * exp(-(v**2)/(b**2))*(1.-f)
			tau_list.append(tau_eff)
			Tb = (Tex * (1.-exp(-tau_eff))) + Tbb * exp(-tau_eff)
			Tb_list.append(Tb)
		def gauss(x, p): # p[0]==mean, p[1]==stdev
			return 1.0/(p[1]*np.sqrt(2*np.pi))*np.exp(-(x-p[0])**2/(2*p[1]**2))	
		p0 = [0,1] # Inital guess is a normal distribution
		errfunc = lambda p, x, y: gauss(x, p) - y # Distance to the target function
		p1, success = opt.leastsq(errfunc, p0[:], args=(vbar_list, Tb_list))
			
		fit_mu, fit_stdev = p1
			
		FWHM = 2*np.sqrt(2*np.log(2))*fit_stdev
		equiv_width = 0
		for Tb in Tb_list:
			equiv_width+=delt_v*(Tb/Tbb-1.)
		plt.text(.880,1.07,r'$\mathbf{^{13}CO\,1\rightarrow0}$ Line Profile', transform=ax.transAxes,horizontalalignment='center',size='18')
		plt.xlabel('Velocity (km/s)',fontweight='bold',fontsize=15)
		plt.ylabel('$\mathbf{T_B(K)}$',fontsize=17)
		plt.text(.900,.090,'Equivalent Width = %4.3f km/s' % equiv_width, transform=ax.transAxes,horizontalalignment='center',weight='bold',size='14')
		#plt.text(.900,.090,'FWHM = %4.3f' % FWHM, transform=ax.transAxes,horizontalalignment='center',weight='bold',size='14')
		plt.plot(vbar_list,Tb_list)
		#plt.plot(vbar_list,tau_list)
		#hist,bins = histogram(clump_vs,bins=50)
		#width = 0.7*(bins[1]-bins[0])
		#center = (bins[:-1]+bins[1:])/2
		#plt.bar(center, hist, align = 'center', width = width)
		plt.draw()


def redraw_Rcl(val):
	Rcl = sRcl.val
	inst.Rcl_recalc(Rcl)
	inst.variables()
	inst.distribute()

def redraw_Rb(val):
	Rb = sRb.val
	inst.Rb_recalc(Rb)
	inst.variables()
	inst.distribute()

def redraw_R(val):
	R = sR.val
	inst.R_recalc(R)
	inst.variables()
	inst.distribute()

def redraw_sigcl(val):
	sigcl = ssigcl.val
	inst.sigcl_recalc(sigcl)
	inst.variables()
	inst.distribute()

def redraw_Vinit(val):
	Vinit = sVinit.val
	inst.Vinit_recalc(Vinit)
	inst.variables()
	inst.distribute()

def redraw_f(val):
	f = sf.val
	inst.f_recalc(f)
	inst.variables()
	inst.distribute()

def redraw_n(val):
	n = sn.val
	inst.n_recalc(n)
	inst.variables()
	inst.distribute()

def redraw_delt_cs(val):
	delt_cs = sdelt_cs.val
	inst.delt_cs_recalc(delt_cs)
	inst.variables()
	inst.distribute()

def redraw_Gamma(val):
	Gamma = sGamma.val
	inst.Gamma_recalc(Gamma)
	inst.variables()
	inst.distribute()

def redraw_gamma(val):
	gamma = sgamma.val
	inst.gamma_recalc(gamma)
	inst.variables()
	inst.distribute()

def redraw_nu(val):
	nu = snu.val
	inst.nu_recalc(nu)
	inst.variables()
	inst.distribute()

def rand_redraw(event):
	inst.variables()
	inst.distribute()

def line_prof_redraw(event):
	inst.variables()
	inst.line_prof()

inst = main()

axN1 = plt.axes([0.10,0.215, 0.45, 0.025])
R_0 = 10
sR = Slider(axN1,'$\mathbf{R_{cloud}}$ (pc)',1,100,valfmt='%2.1f',valinit=R_0)
sR.on_changed(redraw_R)

axN2 = plt.axes([0.10,0.175, 0.45, 0.025])
Rb_0 = 0.05
sRb = Slider(axN2,'$\mathbf{R_{beam}}$ (pc)',0.01,0.1,valfmt='%3.3f',valinit=Rb_0)
sRb.on_changed(redraw_Rb)

axN3 = plt.axes([0.10,0.135, 0.45, 0.025])
Rcl0=.01
sRcl = Slider(axN3,'$\mathbf{R_{cl}}$ (pc)',.005,.1,valfmt='%3.3f',valinit=Rcl0)
sRcl.on_changed(redraw_Rcl)

axN4 = plt.axes([0.10,0.095, 0.45, 0.025])
f_0 = 0.05
sf = Slider(axN4,'f',0.001,1.00,valfmt='%5.5f',valinit=f_0)
sf.on_changed(redraw_f)

axN5 = plt.axes([0.10,0.055, 0.45, 0.025])
sigcl_0 = 1.0
ssigcl = Slider(axN5,'$\mathbf{\sigma_{cl}}$ (km/s)',.1,10,valfmt='%2.1f',valinit=sigcl_0)
ssigcl.on_changed(redraw_sigcl)

axN6 = plt.axes([0.10,0.015, 0.45, 0.025])
Gamma_0 = 1.0
sGamma = Slider(axN6,'$\mathbf{\Gamma_{cl}}$',0.5,10,valfmt='%2.1f',valinit=Gamma_0)
sGamma.on_changed(redraw_Gamma)

inst.variables()
inst.distribute()

axprof = plt.axes([0.71,0.15,0.25,0.080])
plot_line_prof = Button(axprof,'Calculate Line Profile')
plot_line_prof.on_clicked(line_prof_redraw)


axrand = plt.axes([0.71, 0.04, 0.25, 0.080])
randomize = Button(axrand,'Randomize Velocities')
randomize.on_clicked(rand_redraw)

plt.show()
