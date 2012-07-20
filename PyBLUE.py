#!/usr/bin/env python

import sys, os
from numpy import *
from ConfigParser import SafeConfigParser

parser = SafeConfigParser()
parser.read( sys.argv[1] )

#config_general = parser.options( 'general' )
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

unc = []

measurements = array( Nmeasurements * [ 0. ] )
c = 0
for channel in measurements_descriptions:
    measurements[c] = float( parser.get( 'measurements', channel ) )
    c += 1

    unc_ch = {}
    all_unc_this_ch = [ float(num) for num in parser.get( "uncertainties", channel ).split() ]
    u = 0
    for s_unc in uncertainties_descriptions:
        unc_ch[s_unc] = all_unc_this_ch[u]
        u += 1
    unc.append( unc_ch )
print "INFO: Measurements:", measurements


correlations = {}
for s_unc in uncertainties_descriptions:
    m = eval( parser.get( "correlations", s_unc ) )
    correlations[s_unc] = matrix( m )


covtot = matrix( [ [0,0], [0,0] ] )
for syst, rho in correlations.iteritems():
    #if syst == "lumi": continue
    #sv = array( [ unc[0][syst], unc[1][syst] ] )
    slist = []
    for channel in unc:
        slist.append( channel[syst] )
    
    sv = array( slist )
    sm = diag( sv )

    cov = sm * rho * sm

    #cov = matrix( [ [ s[0] * s[0] * rho[0,0], s[0] * s[1] * rho[0,1] ], \
    #      [ s[1] * s[0] * rho[1,0], s[1] * s[1] * rho[1,1] ] ] )
    
    print syst + ":"
    print cov
    covtot = covtot + cov

print "Total Covariance matrix:"
print covtot

invcov = linalg.inv( covtot )

#print "Inverse Total Covariance matrix:"
#print invcov

unitvector = array( Nmeasurements * [ 1 ] )
l = dot( unitvector, invcov )
d = dot( unitvector, invcov )
d = dot( d, unitvector.transpose() )
l = l / d

print "Combination weights:"
print l

res = dot( l, measurements )

totunc = sqrt( dot( l, dot( covtot, l.transpose() ) ) )

print "Result:"
print res, " \pm ", totunc


