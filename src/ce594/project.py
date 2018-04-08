import matplotlib.pyplot as plt
from fe import solve_fe
from math import exp

### Constants
A = 2.0
k = 250.0
l = 10.0
alpha = -5.0


### Mesh
num_elements = 1
num_element_nodes = 6
num_nodes = num_elements * ( num_element_nodes - 1 ) + 1

### Boundary conditions
# Essential
bc_essential = dict()
bc_essential[ num_nodes - 1 ] = 3.0

# Natural
bc_natural = dict()
bc_natural[ 0 ] = -0.2

### Solve
d, x = solve_fe( A, k, l, alpha, num_elements, num_element_nodes, bc_essential, bc_natural )

### Exact solution as a function
def exact ( _xe ):
    return exp( -_xe / 10.0 ) * ( 2.73368 + 0.733676 * exp( _xe / 5.0 ) )

xe = [ float(i) * ( float(l) / ( 250 - 1 ) ) for i in range( 250 ) ]
ye = [ exact( _x ) for _x in xe ]

### Print values we're interested in
print 'Hydraulic head at left end:', d[0]
print 'Flow rate at right end:', ( d[len(d)-1] - d[len(d)-2 ] ) / ( x[len(x)-1] - x[len(x)-2] ) * -k * A

### Plot
plt.plot( x, d, label='FE Solution' )
plt.plot( xe, ye, label='Exact Solution')
plt.legend()
plt.show()