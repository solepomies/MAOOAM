-- aotensor.lua
-- (C) 2013 Lesley De Cruz
-- See LICENSE.txt for license information.

------------------------------------------------------------------------
-- The equation tensor for MAOOAM (the modular arbitrary-order ocean-atmosphere
-- model) which allows for an extensible set of modes in the ocean and in the
-- atmosphere.
------------------------------------------------------------------------

-- Model parameters and the inner products between the truncated fields.
local m = require("params").m
local inprod = require("inprod_analytic")

-- Translate field names into effective coordinates
local natm, noc = inprod.natm, inprod.noc
local function psi(i) return i end
local function theta(i) return i+natm end
local function A(i) return i+2*natm end
local function T(i) return i+2*natm+noc end

--- Kronecker delta
local function kdelta(i,j) return i==j and 1 or 0 end

--- Model tensor. The tensor (in coordinate list form) that holds the
-- coefficients for the bilinear system of differential equations.
-- Constructed as an upper triangular tensor (in the last two dimensions).
local coo_list = {}

--- Coefficient insertion function.
-- Creates an upper triangular tensor (in the last two dimensions)
-- i.e. {i,k,j,v} will be stored as {i,j,k,v} if j<=k.
local function coeff(i,j,k,v)
  if v and v~=0 then
    coo_list[#coo_list+1]= j<=k and {i,j,k,v} or {i,k,j,v}
  end
end

-- Atmosphere equations
------------------------------------------------------------------------

local atmos = inprod.atmos
local a, b, c, d, g, s = atmos.a, atmos.b, atmos.c, atmos.d, atmos.g, atmos.s
local kd, betp = m.kd, m.betp
local kdp, sig0 = m.kdp, m.sig0
local nuap,nuop = m.nuap or 0, m.nuop or 0
local Cpa, Lpa, SBpa, SBpo, sc = m.Cpa, m.Lpa, m.SBpa, m.SBpo, m.sc

coeff(theta(1),0,0, (Cpa/(1 - a[1][1]*sig0))) -- constant forcing
-- Contributions proportional to...
for i=1,natm do
  for j=1,natm do
    -- psi_j
    coeff(psi(i),psi(j),0,-((c[i][j]*betp)/a[i][i]) - (kd*kdelta(i, j))/2 + a[i][j]*nuap)
    coeff(theta(i),psi(j),0, (a[i][j]*kd*sig0)/(-2 + 2*a[i][i]*sig0))
    -- theta_j
    coeff(psi(i),theta(j),0, (kd*kdelta(i, j))/2)
    coeff(theta(i),theta(j),0, (-(sig0*(2*c[i][j]*betp + a[i][j]*(kd + 4*kdp))) + 2*(SBpa + sc*Lpa)*kdelta(i,j))/(-2 + 2*a[i][i]*sig0))
    for k=1,natm do
      -- psi_j x psi_k
      coeff(psi(i),psi(j),psi(k), -(b[i][j][k]/a[i][i]))
      -- theta_j x theta_k
      coeff(psi(i),theta(j),theta(k), -(b[i][j][k]/a[i][i]))
      -- psi_j x theta_k
      coeff(theta(i),psi(j),theta(k), (g[i][j][k] - b[i][j][k]*sig0)/(-1 + a[i][i]*sig0))
      -- theta_j x psi_k
      coeff(theta(i),theta(j),psi(k), (b[i][j][k]*sig0)/(1 - a[i][i]*sig0))
    end
  end
  for j=1,noc do -- A_j
    coeff(psi(i),A(j),0,  kd*d[i][j]/(2*a[i][i]))
    coeff(theta(i),A(j),0, kd*(d[i][j]*sig0)/(2 - 2*a[i][i]*sig0))
    -- T_j
    coeff(theta(i),T(j),0, s[i][j]*(2*SBpo + Lpa)/(2 - 2*a[i][i]*sig0))
  end
end

-- Ocean equations
------------------------------------------------------------------------

local ocean = inprod.ocean
local C, K, M, N, O, W = ocean.C, ocean.K, ocean.M, ocean.N, ocean.O, ocean.W
local dp, rp, G = m.dp, m.rp, m.G
local Cpo, Lpo, sBpo, sBpa = m.Cpo, m.Lpo, m.sBpo, m.sBpa

-- Ocean stream function
-- Contributions propoportional to...
for i=1,noc do
  for j=1,natm do  -- psi_j - theta_j == psi3_j
    coeff(A(i),psi(j), 0, K[i][j]*dp / (M[i][i]+G))
    coeff(A(i),theta(j), 0, -K[i][j]*dp / (M[i][i]+G))
  end
  for j=1,noc do   -- A_j
    coeff(A(i),A(j),   0, -(N[i][j]*betp + M[i][i]*(rp+dp)*kdelta(i,j) - M[i][j]^2*nuop)/(M[i][i]+G))
    for k=1,noc do -- A_j x A_k
      coeff(A(i),A(j),A(k), -C[i][j][k]/(M[i][i]+G))
    end
  end
end

-- Temperature of the ocean
for i=1,noc do
  coeff(T(i),0,0, Cpo*W[i][1])
  for j=1,natm do
    coeff(T(i),theta(j),0, W[i][j]*(2*sc*Lpo + sBpa))
  end
  for j=1,noc do
    coeff(T(i),T(j),0, -(Lpo + sBpo)*kdelta(i,j))
    for k=1,noc do
      coeff(T(i),A(j),T(k), -O[i][j][k])
    end
  end
end

-- @export
return coo_list
