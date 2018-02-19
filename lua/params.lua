-- params.lua
-- (C) 2013-2016 Lesley De Cruz & Jonathan Demaeyer
-- See LICENSE.txt for license information.

------------------------------------------------------------------------
-- Integration parameters and physical parameters for the modular
-- arbitrary-order ocean-atmosphere model MAOOAM.
------------------------------------------------------------------------

local sqrt, cos, sin, pi = math.sqrt, math.cos, math.sin, math.pi

--- Integration parameters
-- @table params.i
local i = {
  t_trans = 1e4,  -- transient period (e.g. 1e7)
  t_run   = 1e5,  -- length of trajectory on the attractor (e.g. 5e8)
  dt      = 1e-2, -- the time step
  out_prefix = "output_maooam_", -- prefix for the output file name
  writeout = 10, -- write out all variables every 'writeout' time units
  statistics = 10, -- accumulate statistics every 'statistics' time units
  compression = true, -- compress output of trajectory (if writeout == true)
  snapshot = 1e4, -- write out state every 'snapshot' time units (false = no writeout)
  walltime = 3600*48, -- run a timer for 'walltime' seconds and exit cleanly when reached.
  -- randomseed = 12345, -- seed for the random generator (integer)
  integrator = "rk2", -- numerical integration scheme, e.g. "rk2" or "rk4"
}

--- Model parameters
-- @table params.m 
local m = {
  --- Scale parameters for the ocean and the atmosphere
  scale = 5e6,     -- the characteristic space scale, L*pi
  f0    = 1.032e-4, -- Coriolis parameter at 45 degrees latitude
  n     = 1.5,     -- aspect ratio (n = 2Ly/Lx ; Lx = 2*pi*L/n; Ly = pi*L)
  rra   = 6370e3,  -- earth radius
  phi0  = pi/4,    -- latitude
  --- Parameters for the ocean
  gp    = 3.1e-2,  -- reduced gravity
  r     = 1e-8,    -- frictional coefficient at the bottom of the ocean
  H     = 5e2,     -- depth of the water layer
  alpha = 0,       -- coefficient characterizing the intensification of the flow on the western boundaries
  d     = 1e-8,    -- the coupling parameter (should be divided by f0 in order to be adimensional)
  nuo   = 0,  -- dissipation in the ocean
  -- Parameters for the atmosphere
  k    = 2e-2, -- bottom friction coefficient
  kp   = 4e-2, -- internal friction coefficient
  sig0 = 1e-1, -- static stability
  nua  = 0,  -- dissipation in the atmosphere
  -- Temperature-related parameters for the ocean
  Go     = 2e8,  -- Specific heat capacity of the ocean (50m layer)
  Co     = 350,  -- Constant short-wave radiation of the ocean
  To0    = 285,-- Stationary solution for the 0-th order ocean temperature
  -- Temperature-related parameters for the atmosphere
  Ga     = 1e7,  -- Specific heat capacity of the atmosphere
  Ca     = 100,  -- Constant short-wave radiation of the atmosphere
  epsa   = 0.76, -- Emissivity coefficient for the grey-body atmosphere
  Ta0    = 270,-- Stationary solution for the 0-th order atmospheric temperature
  -- Other temperature-related parameters/constants
  sc     = 1,    -- Ratio of surface to atmosphere temperature
  lambda = 20,   -- Sensible + turbulent heat exchange between ocean and atmosphere
  R      = 287,  -- Gas constant of dry air
  sB     = 5.6e-8, -- Stefan–Boltzmann constant
}

-- Derived/nondimensionalized quantities
m.L    = m.scale/pi -- characteristic length scale divided by pi
m.LR   = sqrt(m.gp*m.H)/m.f0 -- reduced Rossby deformation radius
m.G    = -m.L^2/m.LR^2    -- \gamma
m.betp = m.L/m.rra*cos(m.phi0)/sin(m.phi0) -- \beta prime
m.rp   = m.r/m.f0         -- r \prime
m.dp   = m.d/m.f0         -- \delta \prime
m.kd   = m.k*2            -- bottom friction coefficient
m.kdp  = m.kp             -- internal friction coefficient
m.nuop = m.nuo/(m.L^2*m.f0) -- dissipation in the ocean
m.nuap = m.nua/(m.L^2*m.f0) -- dissipation in the atmosphere
m.timeunit = 1./(m.f0*24.*3600.) -- dimensional time unit

m.Cpo = m.Co / (m.Go*m.f0) * m.R/(m.f0^2*m.L^2)
m.Lpo = m.lambda / (m.Go*m.f0)
m.Cpa = m.Ca / (m.Ga*m.f0) * m.R/(m.f0^2*m.L^2) / 2 -- Cpa acts on psi1-psi3, not on theta.
m.Lpa = m.lambda / (m.Ga*m.f0)
m.sBpo = 4*m.sB*m.To0^3 / (m.Go*m.f0) -- long wave radiation lost by ocean to atmosphere & space
m.sBpa = 8*m.epsa*m.sB*m.Ta0^3 / (m.Go*m.f0) -- long wave radiation from atmosphere absorbed by ocean
m.SBpo = 2*m.epsa*m.sB*m.To0^3 / (m.Ga*m.f0) -- long wave radiation from ocean absorbed by atmosphere
m.SBpa = 8*m.epsa*m.sB*m.Ta0^3 / (m.Ga*m.f0) -- long wave radiation lost by atmosphere to space & ocean

-- set random seed
require("rand").setrandomseed(i.randomseed, true) -- verbose

--- Generator for output file name
-- Adapt this function if you want to have a custom file name (e.g. if you want
-- it to contain certain parameters)
-- @param suffix Extension, e.g. '.txt' or '.gz'
-- @return output name
function i.getoutfn(suffix)
  return string.format("%sCo%.2e_Ca%.2e_dp%.4e%s.txt",i.out_prefix,m.Co,m.Ca,m.dp,suffix or "")
end

return {i=i, m=m}

