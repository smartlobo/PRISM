# This code generates a PRISM particle from a cloud of vertices. 
# Material properties may be assigned along with orientation, location, and initial velocities (translational and rotational)
# There are also several visualization option (i.e. color and transparency)

# libraries
from yade import pack, geom, qt
from yade.gridpfacet import *
from pylab import *
import numpy as np
import scipy as sp
from scipy import spatial as sp_spatial
from scipy.spatial import Delaunay
import scipy.linalg as la
from scipy.spatial.transform import Rotation as Rt


# gravity
g = 9.81

# enter coordinates for vertices here
#vertices= array([
#			[0.0, 0.0, 0.1],
#			[1.0, 0.0, 0.0],
#			[1.0, 1.0, 0.25],
#			[0.0, 1.0, 0.0],	
#			[0.0, 0.0, 0.80],
#			[1.0, 1.0, 1.0]
					
#					]) 


vertices = array([
			[0.0, 0.0, 0.0],
			[0.0, 1.0, 0.0],
			[1.0, 1.0, 0.0],
			[1.0, 0.0, 0.0],
			[0, 0, 1],
			[1, 0, 1],
			[1, 1, 1],
			[0, 1, 1]
			
			])
					

# material properties
E = 1.0e8 # young's modulus 
phi = 40.0 # angle of friction
rho = 2850.0 # density
R = 0.01 # radius of grid node
nc = 1.0e20 # normal cohesion between grid nodes, grid connections and pfacets
sc = 1.0e20 # shear cohesion between grid nodes, grid connections and pfacets	
etaRoll = 1000 # bending strength
cmatName = 'cmatgravel' # cohesive material name
fmatName = 'fmatgravel' # frictional material name		
		
# change orientation of particle through angles of rotation:
alpha = 0.0 # about z-axis
beta = 0.0 # about y-axis
gamma = 0.0 # about x-axis
		
# velocity components
vx = 0.0
vy = 0.0
vz = 0.0	

# angular velocity vector
omega = Vector3(0.0, 0.0, 0.0)

# colors of PRISM particles
color_pf = [0.4,0.3,0.3]
color_gn = [0.4,0.3,0.3]
color_gc = [0.4,0.3,0.3]

# transparancy and numbering options for visualization
skeleton = False
count_gn = False
skeleton_gc = False
count_gc = False
p_skeleton = False


# function for finding centroid of particle based on original coordinates
def centroid(vertices, R, phi, E, rho, nc, sc, etaRoll, cmatName, fmatName, color_pf, color_gn, color_gc, skeleton, count_gn, skeleton_gc, count_gc, p_skeleton):

	gc_list = []


	O.materials.append( CohFrictMat( young=E,poisson=0.3,density=0.0001, normalCohesion=nc ,shearCohesion=sc, fragile = False, momentRotationLaw=True, isCohesive = True, etaRoll=etaRoll, frictionAngle=phi*pi/180,label=cmatName ) ) 
	
	
	# convex triangulation
	hull = sp_spatial.ConvexHull(vertices)

	d_index = hull.vertices # indices for vertices used in delaunay triangulation

	#print('d_index = ', d_index)

	dvertices = []
	for i in d_index:
		dvertices += [list(vertices[i])]
		
	
	gn_list = []
	for v in dvertices:
		gn_list += [O.bodies.append( gridNode(v,R,wire=skeleton,fixed=False,material=cmatName,color=color_gn, highlight = count_gn) )]

	# convert vertices to an array
	dvertices = np.array(dvertices)

	# convex triangulation using delaunay vertices
	hull = sp_spatial.ConvexHull(dvertices)
	tri = hull.simplices
	
	O.materials.append( FrictMat( young=E,poisson=0.3,density=0.0001,frictionAngle=phi*pi/180,label=fmatName ) )
	
	# create triangulated grid connections 
	for t in tri:
		gna = O.bodies[gn_list[t[0]]]
		gnb = O.bodies[gn_list[t[1]]]
		gnc = O.bodies[gn_list[t[2]]]

		# connect a to b
		if O.interactions.has(gna.id, gnb.id) == False:
				gc_list += [O.bodies.append( gridConnection(gna.id,gnb.id,R,color=color_gc,material=fmatName, wire=skeleton_gc, highlight = count_gc) )]
		# connect b to c	
		if O.interactions.has(gnb.id, gnc.id) == False:
				gc_list += [O.bodies.append( gridConnection(gnb.id,gnc.id,R,color=color_gc,material=fmatName, wire=skeleton_gc, highlight = count_gc) )]
		# connect c to a	
		if O.interactions.has(gnc.id, gna.id) == False:
				gc_list += [O.bodies.append( gridConnection(gnc.id,gna.id,R,color=color_gc,material=fmatName, wire=skeleton_gc, highlight = count_gc) )]


	# Find Delaunay tetrahedron simplices (used for determining centroid of volume) 
	d_pnts = Delaunay(dvertices) # object with attributes
	tetra = d_pnts.simplices # look at simplex attribute


	V_list = [] # volumes of tetrahedra

	# weighted centroids of tetrahedra
	xbarW_list = []
	ybarW_list = []
	zbarW_list = []

	# lists of coordinates of centroids of tetrahedrons
	x_c_list = []
	y_c_list = []
	z_c_list = []

	# get positions of grid nodes of different tetrahedra: a, b, c, and d are vectors which locate the vertices of the tetrahedra in space
	for tet in tetra:
		a = O.bodies[gn_list[int(tet[0])]].state.pos
		b = O.bodies[gn_list[int(tet[1])]].state.pos
		c = O.bodies[gn_list[int(tet[2])]].state.pos
		d = O.bodies[gn_list[int(tet[3])]].state.pos


		# centroids of tetrahedrons
		x_c = (1.0/4.0)*(a[0] + b[0] + c[0] + d[0])
		y_c = (1.0/4.0)*(a[1] + b[1] + c[1] + d[1])
		z_c = (1.0/4.0)*(a[2] + b[2] + c[2] + d[2])
	
		# list of centroid coordinates
		x_c_list += [x_c]
		y_c_list += [y_c]
		z_c_list += [z_c]
		
		
		# calculate volume of i-th tetrahedra
		ad = a - d
		bd = b - d
		cd = c - d

		cross_i = cross(bd, cd)
		dot_i = dot(ad, cross_i)
		
		V_i = abs(dot_i)/6.0

		V_list += [V_i]	

		# calculate weighted centroid of i-th tetrahedra

		xbarW_list += [x_c*V_i]
		ybarW_list += [y_c*V_i]
		zbarW_list += [z_c*V_i]


	# centroid of entire solid polyhedron

	Vol_POLY = sum(V_list) # total volume of polyhedron

	xbar = sum(xbarW_list)/Vol_POLY
	ybar = sum(ybarW_list)/Vol_POLY
	zbar = sum(zbarW_list)/Vol_POLY
	
	bods = gn_list + gc_list
	
	for b in bods:
		O.bodies.erase(b)
		
	return xbar, ybar, zbar


# find centroid of polyhedron defined by vertices
Xbar = centroid(vertices, R, phi, E, rho, nc, sc, etaRoll, cmatName, fmatName, color_pf, color_gn, color_gc, skeleton, count_gn, skeleton_gc, count_gc, p_skeleton)

xbar = Xbar[0]
ybar = Xbar[1]
zbar = Xbar[2]
	
# assign location for center of mass

# to leave the center of mass at its original location, choose xt = Xbar[0], yt = Xbar[1], zt = Xbar[2]
# if you want the center of mass at the origina, choose xt = 0, yt = 0, zt = 0
# or you can choose some other point
xt = Xbar[0]
yt = Xbar[1]
zt = Xbar[2]



# function for rotating and translating vertices
def TrnsRot(vertices, xt, yt, zt, alpha, beta, gamma, xbar, ybar, zbar):				

	# rotations about center of mass
	alpha = alpha*pi/180 # convert to radians
	beta = beta*pi/180
	gamma = gamma*pi/180


	# translate vertices to point of rotation
	vertices = vertices - array( [xbar, ybar, zbar] ) 
		
	# rotate vertex points

	rot_vertices = []


	Rot = array([[ cos(alpha)*cos(beta), cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma) , cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma)], 
			[ sin(alpha)*cos(beta), sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma) , sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma)], 				[ -sin(beta), 		cos(beta)*sin(gamma) , 					  cos(beta)*cos(gamma)]
			])

	for v in vertices:
		rot_vertices += [Rot.dot(v)]

	# redefine vertices
	vertices = rot_vertices

	# translations
	vertices += array( [xt, yt, zt] )
	
	return vertices


# engines for Yade constituitive laws and collision detection
O.engines=[


	ForceResetter(),
	InsertionSortCollider(

		[
		Bo1_Sphere_Aabb(),
		Bo1_Facet_Aabb(),
		Bo1_PFacet_Aabb(),
		Bo1_GridConnection_Aabb(),
		Bo1_Wall_Aabb()
		], 
		

		),
	InteractionLoop(

		[
		Ig2_Sphere_Sphere_ScGeom(),
		Ig2_Facet_Sphere_ScGeom(),
		Ig2_GridConnection_PFacet_ScGeom(),
		Ig2_Sphere_PFacet_ScGridCoGeom(),
		Ig2_GridNode_GridNode_GridNodeGeom6D(),
		Ig2_Sphere_GridConnection_ScGridCoGeom(),
		Ig2_GridConnection_GridConnection_GridCoGridCoGeom(),
		Ig2_PFacet_PFacet_ScGeom()
		],


		[
		Ip2_FrictMat_FrictMat_FrictPhys(),
		Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionNow=True,setCohesionOnNewContacts=False),
		Ip2_FrictViscoMat_FrictViscoMat_FrictViscoPhys(),
		Ip2_FrictMat_FrictViscoMat_FrictViscoPhys()
		],
		

		[
		Law2_ScGeom_FrictPhys_CundallStrack(),
		Law2_ScGeom6D_CohFrictPhys_CohesionMoment(useIncrementalForm = True),
		Law2_ScGridCoGeom_FrictPhys_CundallStrack(), 
		Law2_ScGridCoGeom_CohFrictPhys_CundallStrack(), 
		Law2_GridCoGridCoGeom_FrictPhys_CundallStrack(),
		]),
			
		
		NewtonIntegrator(gravity=(0,0,-g), exactAsphericalRot = True, damping=0.0)
		
		]


# rotate and translate PRISM particle
vertices = TrnsRot(vertices, xt, yt, zt, alpha, beta, gamma, xbar, ybar, zbar)

# calculate products and moments of inertia
def Inertia_tensor(x_c, y_c, z_c, rho, V, a, b, c, d):
	
	x1 = a[0]-x_c; x2 = b[0]-x_c; x3 = c[0]-x_c; x4 = d[0]-x_c 
	y1 = a[1]-y_c; y2 = b[1]-y_c; y3 = c[1]-y_c; y4 = d[1]-y_c	
	z1 = a[2]-z_c; z2 = b[2]-z_c; z3 = c[2]-z_c; z4 = d[2]-z_c

############################################################################### Ixx ##################################################################################

	Ixx = rho*(6.0*V)*(y1**2.0 + y1*y2 + y2**2.0 + y1*y3 + y2*y3 + y3**2.0 + y1*y4 + y2*y4 + y3*y4 + y4**2.0 + z1**2.0 + z1*z2 + z2**2.0 + z1*z3 + z2*z3 + z3**2.0 + z1*z4 + z2*z4 + z3*z4 + z4**2.0)/60.0

############################################################################### Iyy ##################################################################################

	Iyy = rho*(6.0*V)*(x1**2.0 + x1*x2 + x2**2.0 + x1*x3 + x2*x3 + x3**2.0 + x1*x4 + x2*x4 + x3*x4 + x4**2.0 + z1**2.0 + z1*z2 + z2**2.0 + z1*z3 + z2*z3 + z3**2.0 + z1*z4 + z2*z4 + z3*z4 + z4**2.0)/60.0

############################################################################### Izz ##################################################################################

	Izz = rho*(6.0*V)*(x1**2.0 + x1*x2 + x2**2.0 + x1*x3 + x2*x3 + x3**2.0 + x1*x4 + x2*x4 + x3*x4 + x4**2.0 + y1**2.0 + y1*y2 + y2**2.0 + y1*y3 + y2*y3 + y3**2.0 + y1*y4 + y2*y4 + y3*y4 + y4**2.0)/60.0

############################################################################### Iyz ##################################################################################

	Iyz = rho*(6.0*V)*(2.0*y1*z1 + y2*z1 + y3*z1 + y4*z1 + y1*z2 + 2.0*y2*z2 + y3*z2 + y4*z2 + y1*z3 + y2*z3 + 2.0*y3*z3 + y4*z3 + y1*z4 + y2*z4 + y3*z4 + 2.0*y4*z4)/120.0

############################################################################### Ixz ##################################################################################

	Ixz = rho*(6.0*V)*(2.0*x1*z1 + x2*z1 + x3*z1 + x4*z1 + x1*z2 + 2.0*x2*z2 + x3*z2 + x4*z2 + x1*z3 + x2*z3 + 2.0*x3*z3 + x4*z3 + x1*z4 + x2*z4 + x3*z4 + 2.0*x4*z4)/120.0

############################################################################### Ixy ##################################################################################

	Ixy = rho*(6.0*V)*(2.0*x1*y1 + x2*y1 + x3*y1 + x4*y1 + x1*y2 + 2.0*x2*y2 + x3*y2 + x4*y2 + x1*y3 + x2*y3 + 2.0*x3*y3 + x4*y3 + x1*y4 + x2*y4 + x3*y4 + 2.0*x4*y4)/120.0

	return Ixx, Iyy, Izz, Iyz, Ixz, Ixy


# Build PRISM particle
def prism(vertices, R,phi, E, rho, nc, sc, etaRoll, cmatName, fmatName, color_pf, color_gn, color_gc, skeleton, count_gn, skeleton_gc, count_gc, p_skeleton, omega, vx, vy, vz):
	
	gn_list = [] # list of grid nodes
	gc_list = [] # list of grid connections
	pf_list = [] # list of pfacets
		
	ma = 0.001 # added mass for stability
	Ia = 0.001 # added inertia for stability
	
	# grid node material properties
	poisson = 0.3
	alphaKr = 100.0
	alphaKtw = 100.0
	etaRoll = 1.0e8
	etaTwist = 1000
	density = 10.0
	matGN = 'cgn'
	
	O.materials.append( CohFrictMat( young=E,poisson=poisson,density=density, normalCohesion=nc ,shearCohesion=sc, fragile = False, momentRotationLaw=True, isCohesive = True, etaRoll=etaRoll, etaTwist = etaTwist, frictionAngle=phi*pi/180, alphaKr = alphaKr, alphaKtw = alphaKtw,  label=matGN ) ) 
	 
	 
	# grid connection material properties
	poisson = 0.3
	alphaKr = 100.0
	alphaKtw = 100.0
	etaRoll = 1.0e8
	etaTwist = 1000
	nc = 1.0e20
	sc = 1.0e20
	density = 10.0
	matGC = 'fgc'
	
	O.materials.append(FrictMat( young=E,poisson=poisson,density=density, frictionAngle=phi*pi/180, label=matGC ) ) 


	# pfacet material properties
	poisson = 0.3
	alphaKr = 1e100
	alphaKtw = 1e100
	etaRoll = 1.0e8
	etaTwist = 1.0e100
	nc = 1.0e20
	sc = 1.0e20
	density = 10.0
	matPF = 'fpf'
	
	O.materials.append(FrictMat( young=E,poisson=poisson,density=density, frictionAngle=phi*pi/180, label=matPF ) ) 
	
	
	# convex triangulation
	hull = sp_spatial.ConvexHull(vertices)

	d_index = hull.vertices # indices for vertices used in delaunay triangulation


	dvertices = []
	for i in d_index:
		dvertices += [list(vertices[i])]
		

	# create grid nodes from delaunay vertices (dvertices)
	gn_list = []
	for v in dvertices:
		gn_list += [O.bodies.append( gridNode(v,R,wire=skeleton,fixed=False,material=matGN ,color=color_gn, highlight = count_gn) )]

	# convert vertices to an array
	dvertices = np.array(dvertices)

	# convex triangulation using delaunay vertices
	hull = sp_spatial.ConvexHull(dvertices)
	tri = hull.simplices


	# Find Delaunay tetrahedron simplices (used for determining centroid of volume) 
	d_pnts = Delaunay(dvertices) # object with attributes
	tetra = d_pnts.simplices # look at simplex attribute


	V_list = [] # volumes of tetrahedra

	# weighted centroids of tetrahedra
	xbarW_list = []
	ybarW_list = []
	zbarW_list = []

	# lists of coordinates of centroids of tetrahedrons
	x_c_list = []
	y_c_list = []
	z_c_list = []

	# get positions of grid nodes of different tetrahedra: a, b, c, and d are vectors which locate the vertices of the tetrahedra in space
	for tet in tetra:
		a = O.bodies[gn_list[int(tet[0])]].state.pos
		b = O.bodies[gn_list[int(tet[1])]].state.pos
		c = O.bodies[gn_list[int(tet[2])]].state.pos
		d = O.bodies[gn_list[int(tet[3])]].state.pos


		# centroids of tetrahedrons
		x_c = (1.0/4.0)*(a[0] + b[0] + c[0] + d[0])
		y_c = (1.0/4.0)*(a[1] + b[1] + c[1] + d[1])
		z_c = (1.0/4.0)*(a[2] + b[2] + c[2] + d[2])
	
		# list of centroid coordinates
		x_c_list += [x_c]
		y_c_list += [y_c]
		z_c_list += [z_c]
		
		
		# calculate volume of i-th tetrahedra
		ad = a - d
		bd = b - d
		cd = c - d

		cross_i = cross(bd, cd)
		dot_i = dot(ad, cross_i)
		
		V_i = abs(dot_i)/6.0

		V_list += [V_i]


		# calculate weighted centroid of i-th tetrahedra

		xbarW_list += [x_c*V_i]
		ybarW_list += [y_c*V_i]
		zbarW_list += [z_c*V_i]


	# centroid of entire solid polyhedron

	Vol_POLY = sum(V_list) # total volume of polyhedron

	xbar = sum(xbarW_list)/Vol_POLY
	ybar = sum(ybarW_list)/Vol_POLY
	zbar = sum(zbarW_list)/Vol_POLY
	

####################################### centroid has been found, now rotate body using updated vertices (those left over with grid connections) ########################
	

	# create triangulated grid connections 
	for t in tri:
		gna = O.bodies[gn_list[t[0]]]
		gnb = O.bodies[gn_list[t[1]]]
		gnc = O.bodies[gn_list[t[2]]]

		# connect a to b
		if O.interactions.has(gna.id, gnb.id) == False:
				gc_list += [O.bodies.append( gridConnection(gna.id,gnb.id,R,color=color_gc,material=matGC, wire=skeleton_gc, highlight = count_gc) )]
		# connect b to c	
		if O.interactions.has(gnb.id, gnc.id) == False:
				gc_list += [O.bodies.append( gridConnection(gnb.id,gnc.id,R,color=color_gc,material=matGC, wire=skeleton_gc, highlight = count_gc) )]
		# connect c to a	
		if O.interactions.has(gnc.id, gna.id) == False:
				gc_list += [O.bodies.append( gridConnection(gnc.id,gna.id,R,color=color_gc,material=matGC, wire=skeleton_gc, highlight = count_gc) )]

	# create grid node at centroid 
	gn_cntrd = O.bodies.append( gridNode([xbar, ybar, zbar],R,wire=skeleton,fixed=False,material=matGN,color=color_gn, highlight = count_gn) )


	# make vertex grid connections to centroid
	index_list = range(len(gn_list))
	for i in index_list:
		gc_list += [O.bodies.append( gridConnection(gn_list[i],gn_cntrd,R,color=color_gc,material=matGC, wire=skeleton_gc, highlight = count_gc) )]


	# make pfacets
	for t in tri:
		pf_list += [O.bodies.append( pfacet(gn_list[int(t[0])], gn_list[int(t[1])], gn_list[int(t[2])],wire=p_skeleton,material=matPF,color=color_pf) )]

	gn_list += [gn_cntrd] # append centroid grid node to list of grid nodes

	############################################################## calculate inertial properties #############################################################

	V_list = [] # initiate list of volumes of tetrahedrons
	m_list = [] # initiate list of masses of tetrhadedrons

	xi_list = [] # initiate list for x-component of centroids of tetrahedrons
	yi_list = [] # initiate list for y-component of centroids of tetrahedrons
	zi_list = [] # initiate list for z-component of centroids of tetrahedrons


	Ixx = 0.0 # initiate moment of inertia xx
	Iyy = 0.0 # initiate moment of inertia yy
	Izz = 0.0 # initiate moment of inertia zz
	Iyz = 0.0 # initiate product of inertia yz
	Ixz = 0.0 # initiate product of inertia xz
	Ixy = 0.0 # initiate product of inertia xy


	# get positions of grid nodes of different tetrahedra: a, b, c, and d are vectors which locate the vertices of the tetrahedra in space
	for tet in tetra:
		a = O.bodies[gn_list[int(tet[0])]].state.pos
		b = O.bodies[gn_list[int(tet[1])]].state.pos
		c = O.bodies[gn_list[int(tet[2])]].state.pos
		d = O.bodies[gn_list[int(tet[3])]].state.pos

		# centroids of tetrahedrons
		xi = (1.0/4.0)*(a[0] + b[0] + c[0] + d[0])
		yi = (1.0/4.0)*(a[1] + b[1] + c[1] + d[1])
		zi = (1.0/4.0)*(a[2] + b[2] + c[2] + d[2])
	
		# list of centroid coordinates
		xi_list += [xi]
		yi_list += [yi]
		zi_list += [zi]
			
		# calculate volume of i-th tetrahedra
		ad = a - d
		bd = b - d
		cd = c - d

		cross_i = cross(bd, cd)
		dot_i = dot(ad, cross_i)
		
		Vi = abs(dot_i)/6.0

		mi = rho*Vi

		# calculating moments and products of inertia for each tetrahedron abour their respective centroids
		output = Inertia_tensor(xi, yi, zi, rho, Vi, a, b, c, d)
	
		Ixxci = output[0]
		Iyyci = output[1]
		Izzci = output[2]
		Iyzci = output[3]
		Ixzci = output[4]
		Ixyci = output[5]

		V_list += [Vi]

		# use parallel axis theorem for moment of inertia about centroid of entire polyhedron

		Ixx += Ixxci + mi*( (yi - ybar)**2.0 +  (zi - zbar)**2.0 )
		Iyy += Iyyci + mi*( (xi - xbar)**2.0 +  (zi - zbar)**2.0 )
		Izz += Izzci + mi*( (xi - xbar)**2.0 +  (yi - ybar)**2.0 )
		Iyz += Iyzci + mi*( (yi - ybar)*(zi - zbar) )
		Ixz += Ixzci + mi*( (xi - xbar)*(zi - zbar) )
		Ixy += Ixyci + mi*( (xi - xbar)*(yi - ybar) )


	# moment of inertia tensor about centroid of polyhedron
	I = np.array([	[Ixx, -Ixy, -Ixz],
			[-Ixy, Iyy, -Iyz],
			[-Ixz, -Iyz, Izz]

			])


	# sensitivity
	tol = 1e-10

	# need this because of working precision of python will calculate the wrong eigenvectors
	if abs(I[0][1]/I[0][0]) <= tol and abs(I[0][2]/I[0][0]) <= tol:
		if abs(I[1][0]/I[1][1]) <= tol and abs(I[1][2]/I[1][1]) <= tol:
			if abs(I[2][0]/I[2][2]) <= tol and abs(I[2][1]/I[2][2]) <= tol:
				I[0][1] = I[0][2] = I[1][0] = I[1][2] = I[2][0] = I[2][1] = 0.0


	# eigenvalues (inertia about principal axes) and eigenvectors (principal axes) of inertia tensor  
	eig = la.eig(I)

	eigvals = eig[0]
	eigvecs = eig[1]

	I1 = eigvals[0].real
	I2 = eigvals[1].real
	I3 = eigvals[2].real

	E1 = eigvecs[0]
	E2 = eigvecs[1]
	E3 = eigvecs[2]
	
	# sort values for Yade
	I_list = [I1, I2, I3]
	E_list = [[E1[0], E2[0], E3[0]], [E1[1], E2[1], E3[1]], [E1[2], E2[2], E3[2]]] 

	list1 = I_list
	list2 = E_list
	
	list1, list2 = zip(*sorted(zip(list1, list2)))
	
	I11 = list1[0]
	I22 = list1[1]
	I33 = list1[2]	

	E_col1 = np.array(list2[0])
	E_col2 = np.array(list2[1])
	E_col3 = np.array(list2[2])
	
	
	# diagonal eigenvalue matrix (diagonalized inertia tensor in principal coordinates)
	yadeD = np.array([  [I11, 0.0, 0.0],
	                    [0.0, I22, 0.0],
	                    [0.0, 0.0, I33]
	                
	                 ])



	# matrix of eigenvectors
	yadeV = np.array([  [E_col1[0], E_col2[0], E_col3[0]],
	                    [E_col1[1], E_col2[1], E_col3[1]],
	                    [E_col1[2], E_col2[2], E_col3[2]]
	                
	                 ])
	
	
	# determinant of matrix
	det_yadeV = det(yadeV)
	
	# if det == -1 then switch direction of eigenvectors 
	
	if det_yadeV < 0.0:
		E_col1 *= -1.0
		E_col2 *= -1.0
		E_col3 *= -1.0
		
		yadeV = np.array([  [E_col1[0], E_col2[0], E_col3[0]],
	                    [E_col1[1], E_col2[1], E_col3[1]],
	                    [E_col1[2], E_col2[2], E_col3[2]]
	                
	                 ])
	
	# rotation matrix of eigenvectors
	
	Rot = yadeV.transpose() 
	
	r = Rt.from_matrix(Rot)


	# convert rotation matrix to quaternion
	qua = r.as_quat()

	qx = qua[0]
	qy = qua[1]
	qz = qua[2]
	qw = qua[3]

	# set orientation for yade
	orient = Quaternion(qw, qx, qy, qz)
		

	# find contributions of principal moments of inertia from point masses

	# find coordinates in principal axes of centroid separately
	Xp = Vector3(dot(O.bodies[gn_cntrd].state.pos, E1), dot(O.bodies[gn_cntrd].state.pos, E2), dot(O.bodies[gn_cntrd].state.pos, E3))
	xp1 = Xp[0]
	xp2 = Xp[1]
	xp3 = Xp[2]

	pcord = [] # coordinates in principal axes basis
	bodyids = []
	# inertia and mass contribution from vertices
	I11v = 0.0
	I22v = 0.0
	I33v = 0.0
	mv = 0.0

	# finding minimum mass for smallest principal moment of inertia

	r11sqrd = []
	r22sqrd = []
	r33sqrd = []

	# make particle aspherical	
	O.bodies[gn_cntrd].aspherical = True

	# assign moments of inertia about principal axes
	O.bodies[gn_cntrd].state.inertia[0] = yadeD[0][0]  
	O.bodies[gn_cntrd].state.inertia[1] = yadeD[1][1]  
	O.bodies[gn_cntrd].state.inertia[2] = yadeD[2][2] 

	O.bodies[gn_cntrd].state.ori = orient
	
	# check that orientation is correct
	M = O.bodies[gn_cntrd].state.ori.toRotationMatrix() 


	yadeRot = np.array([ [M[0][0], M[0][1], M[0][2]],
			      [M[1][0], M[1][1], M[1][2]],
			      [M[2][0], M[2][1], M[2][2]]
			      
			      ])
			   

	RotT = Rot.transpose()
	A1 = RotT.dot(yadeD)

	yade_check_I = A1.dot(Rot)


	angularMomentum = yade_check_I.dot(omega)
	
	O.bodies[gn_cntrd].state.angMom[0] = angularMomentum[0]
	O.bodies[gn_cntrd].state.angMom[1] = angularMomentum[1]
	O.bodies[gn_cntrd].state.angMom[2] = angularMomentum[2]


	# list of all bodies of polyhedron
	bods = gn_list + gc_list + pf_list


	for b in bods:
		O.bodies[b].state.inertia[0] += Ia 
		O.bodies[b].state.inertia[1] += Ia
		O.bodies[b].state.inertia[2] += Ia

	
	MPOLY = rho*Vol_POLY # total mass of polyhedron
	
	# check total mass
	MT = 0.0
	
	for b in bods:
		MT += O.bodies[b].state.mass # total mass of vertices, corners, and pfacets

	# asssign corrected mass to centroid

	O.bodies[gn_cntrd].state.mass = MPOLY - MT
			
	
	# set stiffness
	for gc in gc_list:
		id1 = O.bodies[gc].shape.node1.id
		id2 = O.bodies[gc].shape.node2.id
		
		O.interactions[id1, id2].phys.kn = 1.0e8 
		O.interactions[id1, id2].phys.kr = 1.0e8 
		O.interactions[id1, id2].phys.ks = 1.0e8 
		O.interactions[id1, id2].phys.ktw = 1.0e8 
		O.interactions[id1, id2].phys.maxRollPl = 1.0e100
		O.interactions[id1, id2].phys.maxTwistPl = 1.0e100
		O.interactions[id1, id2].phys.shearAdhesion =1.0e100
		O.interactions[id1, id2].phys.normalAdhesion =1.0e100
		O.interactions[id1, id2].phys.momentRotationLaw = True
	

	################ SET ANGULAR VELOCITY ##################
	# set velocities (angular and linear)

	O.bodies[gn_cntrd].state.angVel[0] = omega[0] # angular velocity of centroid about x-axis
	O.bodies[gn_cntrd].state.angVel[1] = omega[1] # angular velocity of centroid about y-axis
	O.bodies[gn_cntrd].state.angVel[2] = omega[2] # angular velocity of centroid about z-axis
	
	CG = O.bodies[gn_cntrd].state.pos # centroid of polyhedron
	
	for b in bods:
		r = O.bodies[b].state.pos - CG # position vector going from centroid to other component of body (i.e. grid node or grid connection)
		va = cross(omega, r) # velocity from angular contribution 
		vl = Vector3(vx, vy, vz) # velocity from linear contribution
		vt = va + vl # total velocity
		
		O.bodies[b].state.vel = vt
	
	return gn_cntrd, MPOLY, I11, I22, I33, I, gn_list, Vol_POLY, E1, E2, E3, xbar, ybar, zbar, yadeD, bods
	
	
prism(vertices, R,phi, E, rho, nc, sc, etaRoll, cmatName, fmatName, color_pf, color_gn, color_gc, skeleton, count_gn, skeleton_gc, count_gc, p_skeleton, omega, vx, vy, vz)
