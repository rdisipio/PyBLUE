#!/usr/bin/env python

import sys, os
from numpy import *
from ConfigParser import SafeConfigParser


############################################


def CalcCovariance( variation = 'm' ):
    covtot = matrix( [ [0] * Nobservables ] * Nobservables )
    for syst, rho in correlations.iteritems():
        #if syst == "lumi": continue
        #sv = array( [ unc[0][syst], unc[1][syst] ] )
        slist = []
        for channel in unc[variation]:
            slist.append( channel[syst] )

        sv = array( slist )
        sm = diag( sv )

        cov = sm * rho * sm

        print "Source of uncertainty:",syst
        print cov
        covtot = covtot + cov

    print "Total Covariance matrix for variation", variation, ":"
    print covtot
    return covtot

############################################

if len( sys.argv ) == 1:
    print "Usage: ./PyBLUE config.blue"
    exit(0)
    
parser = SafeConfigParser()
parser.read( sys.argv[1] )

Nobservables = int( parser.get( 'general', 'observables' ) )
print "INFO: Number of observables:", Nobservables

measurements_descriptions = parser.get( 'general', 'measurements' ).split()
Nmeasurements = len( measurements_descriptions )
print "INFO: Number of channels:", Nmeasurements
print "INFO: Channels:", measurements_descriptions

uncertainties_descriptions = parser.get( 'general', 'uncertainties' ).split()
Nuncertainties = len( uncertainties_descriptions )
print "INFO: Number of uncertainties:", Nuncertainties
print "INFO: Uncertainties:", uncertainties_descriptions

unc  = {
    'u' : [],
    'm' : [],
    'd' : []
    }

blue_central = {}
blue_unc     = {}
for variation in unc.keys():
    blue_central[variation] = 0.
    blue_unc[variation]     = 0.

measurements = array( Nmeasurements * [ 0. ] )
n_c = 0
for channel in measurements_descriptions:
    measurements[n_c] = float( parser.get( 'measurements', channel ) )
    n_c += 1

    unc_ch = {
        'u' : {},
        'm' : {},
        'd' : {}
        }
    channel_u = channel + "_u"
    channel_d = channel + "_d"
    
    all_unc_this_ch_u = [ float(num) for num in parser.get( "uncertainties", channel_u ).split() ]
    all_unc_this_ch_d = [ float(num) for num in parser.get( "uncertainties", channel_d ).split() ]
    #all_unc_this_ch_m = [ 0.5 * ( up + down ) for up in all_unc_this_ch_u for down in all_unc_this_ch_d ]
    all_unc_this_ch_m = []
    for nx in range( len(all_unc_this_ch_u) ):
        all_unc_this_ch_m.append( 0.5 * (all_unc_this_ch_u[nx] + all_unc_this_ch_d[nx])  )

    n_u = 0
    for s_unc in uncertainties_descriptions:
#        print s_unc, n_u, channel, all_unc_this_ch_u[n_u]
        unc_ch['u'][s_unc] = all_unc_this_ch_u[n_u]
        unc_ch['m'][s_unc] = all_unc_this_ch_m[n_u]
        unc_ch['d'][s_unc] = all_unc_this_ch_d[n_u]
        n_u += 1

    for variation in unc_ch.keys():
        unc[variation].append( unc_ch[variation] )

        
print "INFO: Measurements:", measurements
print

correlations = {}
for s_unc in uncertainties_descriptions:
    m = eval( parser.get( "correlations", s_unc ) )
    correlations[s_unc] = matrix( m )


for variation in [ 'u', 'd', 'm' ]:
    print "Variation", variation
    cov = CalcCovariance( variation )

    # pinv = pseudo-inverse. workaround for non-invertible matrices
    #invcov = linalg.inv( cov )
    invcov = linalg.pinv( cov )
    print "Inverted covariance matrix:"
    print invcov

    unitvector = array( Nmeasurements * [ 1 ] )

    l = dot( unitvector, invcov )
    d = dot( unitvector, invcov )
    d = dot( d, unitvector.transpose() )
#    l = dot( unitvector.transpose(), invcov )
#    d = dot( unitvector.transpose(), invcov )
#    d = dot( d, unitvector )

    l = l / d

    print "Combination weights for variation", variation, ":"
    print l

    res = dot( l, measurements )

    totunc = sqrt( dot( l, dot( cov, l.transpose() ) ) )

    blue_central[variation] = float( res )
    blue_unc[variation]     = float( totunc )
    print
    
print
print "/------------------------------------/"
print
print "Final results:"
print
for variation in [ 'u', 'd', 'm' ]:
    print "Combination(%s) = %5.4f \pm %10.7f" % ( variation, blue_central[variation], blue_unc[variation] ) 
    print
    
# AIB
#R_u = blue_unc['u'] / ( blue_unc['u'] + blue_unc['d'] )
#sigma_u = 2. * R_u * blue_unc['m']
#sigma_d = 2. * ( 1. - R_u ) * blue_unc['m']
#
#print "Combination(AIB) = %5.4f +%5.4f -%5.4f" % ( blue_central['m'], sigma_u, sigma_d )

denom = 0.5 * ( blue_unc['u'] + blue_unc['d'] )
factor = blue_unc['m'] / denom

upperError = factor * blue_unc['u']
lowerError = factor * blue_unc['d']

print "Combination(AIB) = %10.7f +%10.7f -%10.7f" % ( blue_central['m'], upperError, lowerError )
#print "Combination(AIB) = %10.7f +%8.7f -%8.7f" % ( blue_central['m'], blue_unc['u'], blue_unc['d'] )
