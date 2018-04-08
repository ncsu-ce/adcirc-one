import matplotlib.pyplot as plt

def shape_functions ( num_element_nodes ):

    # Create xi array
    xi = [ -1. + float(i) * 2. / float( num_element_nodes - 1 )  for i in range( num_element_nodes ) ]

    # Create N and dN arrays
    N = []
    dN = []

    for i in range( num_element_nodes ):

        # Create shape function
        f = []

        for j in range( num_element_nodes ):

            if i != j:

                f.append( lambda _xi, _i=xi[i], _j=xi[j]: ( _xi - _j ) / ( _i - _j ) )

        N.append( lambda _xi, _f=f: reduce( lambda a,b: a*b, map( lambda c: c(_xi), _f ) ) )


        # Create derivative of shape function
        mf = []

        for j in range( num_element_nodes ):

            if i != j:

                m = 1.0 / ( xi[ i ] - xi[ j ] )
                f = []

                for k in range( num_element_nodes ):

                    if k != i and k != j:

                        f.append( lambda _xi, _i=xi[ i ], _k=xi[ k ]: ( _xi - _k ) / ( _i - _k ) )

                    else:

                        f.append( lambda _xi: 1.0 )

                mf.append( ( m, f ) )

        dN.append( lambda _xi, _mf=mf: reduce( lambda a,b: a+b, [ _m * reduce( lambda c,d: c*d, map( lambda __f: __f(_xi), _f ) ) for (_m,_f) in _mf ] ) )

    return N, dN, xi

def plot_shape_functions ( shape_fns ):

    npoints = 250
    x = [ -1 + i * 2. / ( npoints- 1 )  for i in range( npoints ) ]

    for f in shape_fns:
        plt.plot( x, [ f( i ) for i in x ] )

    plt.show()

def plot_f_df ( f, df ):

    npoints = 250
    x = [ -1 + i * 2. / ( npoints- 1 )  for i in range( npoints ) ]

    plt.plot( x, [ f( i ) for i in x ] )
    plt.plot( x, [ df( i ) for i in x ] )
    plt.show()


if __name__ == '__main__':

    nen = 7
    N, dN, xi = shape_functions( nen )

    # Check values
    for i in range( nen ):
        print 'N' + str(i) + '(', xi[i], ') = ', N[i](xi[i])

    plot_shape_functions( dN )
    # plot_f_df( N[0], dN[0] )

