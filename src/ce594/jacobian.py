# Tranform 1d coordinate from normalized value on [-1,1] to value on [xl,xr]
def denormalize ( xi, xl, xr ):
    return 0.5 * ( xr + xl ) + ( xi / 2.0 ) * ( xr - xl )

# The derivative of the denormalize function
def jacobian ( xl, xr ):
    return 0.5 * ( xr - xl )

# Generate a lambda function for a coordinate range for the denormalize function
def fdenormalize ( xl, xr ):
    return lambda xi: denormalize( xi, xl, xr )