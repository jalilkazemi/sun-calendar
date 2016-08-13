##### Sun and shadow's trajectory
library(ggplot2)
SunPos <- function(tilt, alt) {
  DegreeToRadiant <- function(x) x/180*pi
  gma <- cos(DegreeToRadiant(tilt))
  lmda <- cos(DegreeToRadiant(alt))
  
  minTheta <- max(0, gma*lmda - sqrt((1-gma^2)*(1-lmda^2)))
  maxTheta <- max(0, gma*lmda + sqrt((1-gma^2)*(1-lmda^2)))
  theta <- seq(minTheta, maxTheta, length.out = 100)
  tau <- (gma - lmda*theta) / sqrt(1-theta^2) / sqrt(1-lmda^2)
  tau[tau < -1] <- -1
  tau[tau > 1] <- 1
  
  half <- data.frame(theta=pi/2-acos(theta), tau=acos(tau))
  Rev <- function(x) x[length(x):1]
  full <- rbind(half, data.frame(theta=Rev(half$theta), tau=-Rev(half$tau)))
  
  RemoveJump <- function(x1, x2) floor((x1-x2)/(2*pi)+0.5)*2*pi
  discontPoint <- nrow(half)
  shift <- RemoveJump(full$tau[discontPoint], full$tau[discontPoint+1])
  full$tau[discontPoint+seq_len(discontPoint)] <- shift + full$tau[discontPoint+seq_len(discontPoint)]
  
  shadowSize <- cos(full$theta)/sin(full$theta)
  full$x <- sin(pi+full$tau)*shadowSize
  full$y <- cos(pi+full$tau)*shadowSize
  
  RadiantToDegree <- function(x) x/pi*180
  full$theta <- RadiantToDegree(full$theta)
  full$tau <- RadiantToDegree(full$tau)
  
  names(full) <- c('altitude', 'orientation', 'x', 'y')
  finiteView <- !is.infinite(full$x) & !is.infinite(full$y)
  return(list(data=full,
              sun=ggplot(data = full, aes(x = orientation, y = altitude)) + 
                ggtitle(sprintf('The trajectory of Sun on sky dome\n tilt=%.1f, alt=%.1f', tilt, alt)) + 
                geom_path(),
              shadow=ggplot(data = subset(full, finiteView), aes(x = x, y = y)) + 
                ggtitle(sprintf('The trajectory of shadow\n tilt=%.1f, alt=%.1f', tilt, alt)) + 
                geom_path()
              ))
}

minTilt <- 90 - 23.5
maxTilt <- 90 + 23.5

SunPos(90, 40)$sun
SunPos(minTilt, 40)$sun
SunPos(maxTilt, 40)$sun

SunPos(90, 20)$sun
SunPos(minTilt, 20)$sun
SunPos(maxTilt, 20)$sun

SunPos(90, 180-40)$sun
SunPos(minTilt, 180-40)$sun
SunPos(maxTilt, 180-40)$sun

SunPos(90, 180-20)$sun
SunPos(minTilt, 180-20)$sun
SunPos(maxTilt, 180-20)$sun

SunPos(90, 40)$shadow
SunPos(minTilt, 40)$shadow
SunPos(maxTilt, 40)$shadow

SunPos(90, 20)$shadow
SunPos(minTilt, 20)$shadow
SunPos(maxTilt, 20)$shadow

SunPos(90, 180-40)$shadow
SunPos(minTilt, 180-40)$shadow
SunPos(maxTilt, 180-40)$shadow

SunPos(90, 180-20)$shadow
SunPos(minTilt, 180-20)$shadow
SunPos(maxTilt, 180-20)$shadow


#### normal to shadow trajectory
ShadowPosImplicit <- function(tilt, alt) {
  DegreeToRadiant <- function(x) x/180*pi
  gma <- cos(DegreeToRadiant(tilt))
  lmda <- cos(DegreeToRadiant(alt))
  A <- gma / sqrt(1-lmda^2)
  B <- lmda / sqrt(1-lmda^2)

  xs <- seq(-20, 20, by =.01)
  ys <- seq(-20, 20, by =.01)
  x <- matrix(rep(xs, times = length(ys)), nrow = length(xs))
  y <- matrix(rep(ys, each = length(ys)), nrow = length(xs))
  z <- y + A*sqrt(x^2+y^2+1) - B
  contour(x = xs, y = ys, z = z, levels = 0)
}

ShadowPosImplicit(90, 40)
ShadowPosImplicit(90, 20)
ShadowPosImplicit(minTilt, 40)
ShadowPosImplicit(minTilt, 20)
ShadowPosImplicit(minTilt, 180-40)
ShadowPosImplicit(maxTilt, 180-20)
