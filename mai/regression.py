import numpy as np

def OLS_fit(xyz):
	"""
	Make ordinary least squares fit to z=a+bx+cy and return the normal vector

	Args:
		xyz (numpy array): 2D numpy array of XYZ values (N rows, 3 cols)
	Returns:
		normal_vec (numpy array): 1D numpy array for the normal vector
	"""
	#Get fit parameters (a,b,c)
	x = xyz[:,0][np.newaxis].T
	y = xyz[:,1][np.newaxis].T
	z = xyz[:,2][np.newaxis].T
	onevec = np.ones((len(x),1))
	A = np.hstack((x,y,onevec))
	B = z
	fit = np.squeeze(np.dot(np.linalg.pinv(np.dot(A.T,A)),np.dot(A.T,B)))
	
	#Test quality of fit
	z_fit = fit[0]*x+fit[1]*y+fit[2]
	ss_res = sum((z_fit-z)**2)[0]
	ss_tot = sum((z-np.mean(z))**2)[0]
	
	#Calculate r^2
	r2 = 1-ss_res/ss_tot
	normal_vec = np.array([fit[0],fit[1],-1])
	if r2 < (1-10**-14) and len(x) == 2:
		raise ValueError('Poor linear fit to two points?!')

	return normal_vec

def TLS_fit(xyz):
	"""
	Make total least squares fit to ax+by+cz+d=0 and return the normal vector

	Args:
		xyz (numpy array): 2D numpy array of XYZ values (N rows, 3 cols)
	Returns:
		rmse (float): root mean square error of fit
		normal_vec (numpy array): 1D numpy array for the normal vector
	"""
	#Use SVD method to perform TLS regression. TLS is recommended over OLS
	#because OLS will fail to fit a perfectly or nearly vertical plane,
	#whereas TLS will not have this issue

	#Get fit parameters
	xyz_mean = np.mean(xyz,axis=0)
	xyz_sub = xyz-xyz_mean
	[u,s,v] = np.linalg.svd(xyz_sub,full_matrices=False)
	v = v.T
	normal_vec = v[:,-1]
	a = normal_vec[0]
	b = normal_vec[1]
	c = normal_vec[2]
	d = -np.dot(xyz_mean,normal_vec)

	#Test quality of fit
	fit = a*xyz[:,0]+b*xyz[:,1]+c*xyz[:,2]+d
	ss_res = sum(fit**2)
	rmse = (ss_res/np.shape(xyz)[0])**0.5

	return rmse, normal_vec
