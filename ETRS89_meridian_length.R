# calculates the length of an ellipsoidal arc, from the equator to latitude lat_in, 
# along an ellipse of longitude, written in geodetic coordinate of latitude
# lat_in should be in degrees
ell_arc <- function(lat_in, a, b){
  lat <- lat_in * pi / 180
  n <- (a-b)/(a+b)
  arc <- (1 + n^2/4 + n^4/64) * lat + (-3*n/2 + 3*n^3/16) * sin(2*lat) +
         (15*n^2/16 - 15*n^4/64) * sin(4*lat) + (-35*n^3/48) * sin(6*lat) + 
         (315*n^4/512) * sin(8*lat)
  return( (a+b)*arc/2 )
}

# calculates the distance between two points on the ETRS89
# ellipsoid located on the sam longitude, i.e. 
# two points (long, lat_1) and (long, lat_2) along an ellipse of longitude
# lat_1, lat_2 should be in degrees
ell_dist <- function(lat_1, lat_2, a, b){
  return(ell_arc(lat_2, a, b) - ell_arc(lat_1, a, b))
}

# derivative of the distance along an ellipse of longitude
# lat_in is in degrees
der_ell_dist <- function(lat_in, a, b){
  e_2 <- 1 - b^2/a^2
  denom <- 1 - e_2*(sin( lat_in * pi/180 ))^2
  return( a*(1-e_2) / ( denom*sqrt(denom)) )
}

# given a point (long, lat) on the ETRS89 ellipsoid, 
# find the point (long, lat_h) along the ellips of longitude
# such that the elliptic distance between the two points is
# prescribed as h. Found via Newton's method
# lat is in degrees
find_lat <- function(lat, h, a, b){
  lat_h <- lat
  el_dst <- -h
  while(abs(el_dst) > 0.000001){
    lat_h <- lat_h - el_dst / der_ell_dist(lat_h, a, b)
    el_dst <- ell_dist(lat, lat_h, a, b) - h
  }
  return(lat_h)
}

#ETRS89 semi-major and semi-minor axes:
a = 6378137.000
b = 6356752.314140
