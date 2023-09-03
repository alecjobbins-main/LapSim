
% This script generates a spline describing the racetrack
% It then calculates the curvature as a function of arclength (distance)

%% Map Making Instructions
% insert image of track into solidworks part file
% create style spline in Solidworks (degree 5 seems to work well) on track
% export points by running macro Sourbh-ExtractSplinePoints
% select sketch with spline then click tools->macro->run->select file

% open excel file in Matlab and import X and Y coordinates
% if you didn't click the spline points in order on Solidworks, they may be out of order
% to fix the order, scatter(X,Y) and use the brush to select the first point
% right-click 'Copy Data to Clipboard' and  paste it into an excel file
% continue with the rest of the points in order

% import the re-ordered X and Y coordinates as column vectors
% set degree as one greater than you used for drawing the spline 
% make sure to add nurbs_toolbox to search path
% run this file

%% NURBS Toolbox
%load('bezierpoints.mat');
clc;
clear all;
close all;

% spline tools from NURBS Toolbox by D.M. Spink
X = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 15 14 13 12 11 10 9 9 10 11 12 13 14 15 16 17 16 15 14];
Y = [0 1 3 6 10 15 21 28 36 45 55 65 74 83 90 75 60 45 40 45 30 45 30 20 10 5 1 10 10 5 -1 -10 -20 -1 -30];
coefs = [(X); (Y)]; 
degree = 6;
knots = [zeros(1,degree-1) linspace(0,1,numel(X)-degree+2) ones(1,degree-1)];

% creates NURBS (Non-uniform rational basis spline)
nurbs = nrbmak(coefs,knots);


subdivisions = 1002;
p = nrbeval(nurbs,linspace(0.0,1.0,subdivisions)); 
% evaluating curvature
dY = diff(p(2,:))./diff(linspace(0,1,subdivisions));   % first derivative
ddY = diff(dY)./diff(linspace(0,1,subdivisions-1)); %second derivative
dX = diff(p(1,:))./diff(linspace(0,1,subdivisions)); %first derivative
ddX = diff(dX)./diff(linspace(0,1,subdivisions-1)); %second derivative
dX = dX(1:end-1);
dY = dY(1:end-1);

curvature = (dX.*ddY-ddX.*dY)./(dX.*dX+dY.*dY).^(3/2);

% Plotting
figure();
for i = 1:length(p(1,:))-2
    if curvature(i) > 0
        plot(p(1,i),p(2,i), 'g.');
    else
        plot(p(1,i),p(2,i), 'r.');
    end
    hold on;
end
plot(X,Y,'o');
grid on;

% evaluating arclength
seglength = sqrt(sum(diff(p,[],2).^2,1));
total_arclength = sum(seglength);
arclength = linspace(0,total_arclength,numel(curvature));

% plotting
figure(2)
plot(arclength,curvature);
xlabel('Distance','FontSize',15);
ylabel('Curvature','FontSize',15);

% when satisfied with results, save arclength and curvature data
save('track_endurance_2019.mat','arclength','curvature');

function nurbs = nrbmak(coefs,knots) 
% 
% Function Name: 
%  
%   nrbmak - Construct the NURBS structure given the control points 
%            and the knots. 
%  
% Calling Sequence: 
%  
%   nurbs   = nrbmak(cntrl,knots); 
%  
% Parameters: 
%  
%   cntrl       : Control points, these can be either Cartesian or 
% 		homogeneous coordinates. 
%  
% 		For a curve the control points are represented by a 
% 		matrix of size (dim,nu) and for a surface a multidimensional 
% 		array of size (dim,nu,nv). Where nu is number of points along 
% 		the parametric U direction, and nv the number of points 
%               along the V direction. Dim is the dimension valid options 
% 		are 
% 		2 .... (x,y)        2D Cartesian coordinates 
% 		3 .... (x,y,z)      3D Cartesian coordinates 
% 		4 .... (wx,wy,wz,w) 4D homogeneous coordinates 
%  
%   knots	: Non-decreasing knot sequence spanning the interval 
%               [0.0,1.0]. It's assumed that the curves and surfaces 
%               are clamped to the start and end control points by knot 
%               multiplicities equal to the spline order. 
%               For curve knots form a vector and for a surface the knot 
%               are stored by two vectors for U and V in a cell structure. 
%               {uknots vknots} 
%                
%   nurbs 	: Data structure for representing a NURBS curve. 
%  
% NURBS Structure: 
%  
%   Both curves and surfaces are represented by a structure that is 
%   compatible with the Spline Toolbox from Mathworks 
%  
% 	nurbs.form   .... Type name 'B-NURBS' 
% 	nurbs.dim    .... Dimension of the control points 
% 	nurbs.number .... Number of Control points 
%       nurbs.coefs  .... Control Points 
%       nurbs.order  .... Order of the spline 
%       nurbs.knots  .... Knot sequence 
%  
%   Note: the control points are always converted and stored within the 
%   NURBS structure as 4D homogeneous coordinates. A curve is always stored  
%   along the U direction, and the vknots element is an empty matrix. For 
%   a surface the spline degree is a vector [du,dv] containing the degree 
%   along the U and V directions respectively. 
%  
% Description: 
%  
%   This function is used as a convenient means of constructing the NURBS 
%   data structure. Many of the other functions in the toolbox rely on the  
%   NURBS structure been correctly defined as shown above. The nrbmak not 
%   only constructs the proper structure, but also checks for consistency. 
%   The user is still free to build his own structure, in fact a few 
%   functions in the toolbox do this for convenience. 
%  
% Examples: 
%  
%   Construct a 2D line from (0.0,0.0) to (1.5,3.0). 
%   For a straight line a spline of order 2 is required. 
%   Note that the knot sequence has a multiplicity of 2 at the 
%   start (0.0,0.0) and end (1.0 1.0) in order to clamp the ends. 
%  
%   line = nrbmak([0.0 1.5; 0.0 3.0],[0.0 0.0 1.0 1.0]); 
%   nrbplot(line, 2); 
%  
%   Construct a surface in the x-y plane i.e 
%      
%     ^  (0.0,1.0) ------------ (1.0,1.0) 
%     |      |                      | 
%     | V    |                      | 
%     |      |      Surface         | 
%     |      |                      | 
%     |      |                      | 
%     |  (0.0,0.0) ------------ (1.0,0.0) 
%     | 
%     |------------------------------------> 
%                                       U  
% 
%   coefs = cat(3,[0 0; 0 1],[1 1; 0 1]); 
%   knots = {[0 0 1 1]  [0 0 1 1]} 
%   plane = nrbmak(coefs,knots); 
%   nrbplot(plane, [2 2]); 
 
%  D.M. Spink 
%  Copyright (c) 2000. 
 
nurbs.form   = 'B-NURBS'; 
nurbs.dim    = 4; 
np = size(coefs); 
dim = np(1); 
if iscell(knots) 
  % constructing a surface  
  nurbs.number = np(2:3); 
  if (dim < 4) 
    nurbs.coefs = repmat([0.0 0.0 0.0 1.0]',[1 np(2:3)]); 
    nurbs.coefs(1:dim,:,:) = coefs;   
  else 
    nurbs.coefs = coefs; 
  end 
  uorder = size(knots{1},2)-np(2); 
  vorder = size(knots{2},2)-np(3); 
  uknots = sort(knots{1}); 
  vknots = sort(knots{2}); 
  uknots = (uknots-uknots(1))/(uknots(end)-uknots(1)); 
  vknots = (vknots-vknots(1))/(vknots(end)-vknots(1)); 
  nurbs.knots = {uknots vknots}; 
  nurbs.order = [uorder vorder]; 
 
else 
 
  % constructing a curve 
  nurbs.number = np(2); 
  if (dim < 4) 
    nurbs.coefs = repmat([0.0 0.0 0.0 1.0]',[1 np(2)]); 
    nurbs.coefs(1:dim,:) = coefs;   
  else 
    nurbs.coefs = coefs; 
  end 
  nurbs.order = size(knots,2)-np(2); 
  knots = sort(knots); 
  nurbs.knots = (knots-knots(1))/(knots(end)-knots(1)); 
 
end 
end

function [p,w] = nrbeval(nurbs,tt) 
%  
% Function Name: 
%  
%   nrbeval - Evaluate a NURBS at parameteric points 
%  
% Calling Sequence: 
%  
%   [p,w] = nrbeval(crv,ut) 
%   [p,w] = nrbeval(srf,{ut,vt}) 
%  
% Parameters: 
%  
%   crv		: NURBS curve, see nrbmak. 
%  
%   srf		: NURBS surface, see nrbmak. 
%  
%   ut		: Parametric evaluation points along U direction. 
% 
%   vt		: Parametric evaluation points along V direction. 
%  
%   p		: Evaluated points on the NURBS curve or surface as cartesian 
% 		coordinates (x,y,z). If w is included on the lhs argument list 
% 		the points are returned as homogeneous coordinates (wx,wy,wz). 
%  
%   w		: Weights of the homogeneous coordinates of the evaluated 
% 		points. Note inclusion of this argument changes the type  
% 		of coordinates returned in p (see above). 
%  
% Description: 
%  
%   Evaluation of NURBS curves or surfaces at parametric points along the  
%   U and V directions. Either homogeneous coordinates are returned if the  
%   weights are requested in the lhs arguments, or as cartesian coordinates. 
%   This function utilises the 'C' interface bspeval. 
%  
% Examples: 
%  
%   Evaluate the NURBS circle at twenty points from 0.0 to 1.0 
%  
%   nrb = nrbcirc; 
%   ut = linspace(0.0,1.0,20); 
%   p = nrbeval(nrb,ut); 
%  
% See: 
%   
%     bspeval 
% 
 
%  D.M. Spink 
%  Copyright (c) 2000. 
 
if nargin < 2 
  error('Not enough input arguments'); 
end 
 
foption = 1;    % output format 3D cartesian coordinates 
if nargout == 2 
  foption = 0;  % output format 4D homogenous coordinates  
end 
    
if ~isstruct(nurbs) 
  error('NURBS representation is not structure!'); 
end 
 
if ~strcmp(nurbs.form,'B-NURBS') 
  error('Not a recognised NURBS representation'); 
end 
 
if iscell(nurbs.knots) 
  % NURBS structure represents a surface 
 
  num1 = nurbs.number(1); 
  num2 = nurbs.number(2); 
  degree = nurbs.order-1; 
 
  if iscell(tt) 
    % Evaluate over a [u,v] grid 
    % tt{1} represents the u direction 
    % tt{2} represents the v direction 
 
    nt1 = length(tt{1}); 
    nt2 = length(tt{2}); 
     
    % Evaluate along the v direction 
    val = reshape(nurbs.coefs,4*num1,num2); 
    val = bspeval(degree(2),val,nurbs.knots{2},tt{2}); 
    val = reshape(val,[4 num1 nt2]); 
     
    % Evaluate along the u direction 
    val = permute(val,[1 3 2]); 
    val = reshape(val,4*nt2,num1); 
    val = bspeval(degree(1),val,nurbs.knots{1},tt{1}); 
    val = reshape(val,[4 nt2 nt1]); 
    val = permute(val,[1 3 2]); 
 
    w = val(4,:,:); 
    p = val(1:3,:,:); 
    if foption 
      p = p./repmat(w,[3 1 1]); 
    end 
 
  else 
 
    % Evaluate at scattered points 
    % tt(1,:) represents the u direction 
    % tt(2,:) represents the v direction 
 
    nt = size(tt,2); 
 
    val = reshape(nurbs.coefs,4*num1,num2); 
    val = bspeval(degree(2),val,nurbs.knots{2},tt(2,:)); 
    val = reshape(val,[4 num1 nt]); 
 
 
    % evaluate along the u direction 
    pnts = zeros(4,nt); 
    for v = 1:nt 
      coefs = squeeze(val(:,:,v)); 
      pnts(:,v) = bspeval(degree(1),coefs,nurbs.knots{1},tt(1,v)); 
    end 
 
    w = pnts(4,:); 
    p = pnts(1:3,:); 
    if foption 
       p = p./w; 
    end 
         
  end 
 
else 
 
  % NURBS structure represents a curve 
  %  tt represent a vector of parametric points in the u direction 
 
  val = bspeval(nurbs.order-1,nurbs.coefs,nurbs.knots,tt);    
 
  w = val(4,:); 
  p = val(1:3,:); 
  if foption 
    p = p./repmat(w,3,1); 
  end 
 
end 
end


function p = bspeval(d,c,k,u) 
%  
% Function Name: 
%  
%   bspeval - Evaluate a univariate B-Spline. 
%  
% Calling Sequence: 
%  
%   p = bspeval(d,c,k,u) 
%  
% Parameters: 
%  
%   d	: Degree of the B-Spline. 
%  
%   c	: Control Points, matrix of size (dim,nc). 
%  
%   k	: Knot sequence, row vector of size nk. 
%  
%   u	: Parametric evaluation points, row vector of size nu. 
%  
%   p	: Evaluated points, matrix of size (dim,nu) 
%  
% Description: 
%  
%   Evaluate a univariate B-Spline. This function provides an interface to 
%   a toolbox 'C' routine. 
nu = numel(u); 
[mc,nc] = size(c); 
                                                %   int bspeval(int d, double *c, int mc, int nc, double *k, int nk, double *u,int nu, double *p){ 
                                                %   int ierr = 0; 
                                                %   int i, s, tmp1, row, col; 
                                                %   double tmp2; 
                                                % 
                                                %   // Construct the control points 
                                                %   double **ctrl = vec2mat(c,mc,nc); 
                                                % 
                                                %   // Contruct the evaluated points 
p = zeros(mc,nu);                               %   double **pnt = vec2mat(p,mc,nu); 
                                                % 
                                                %   // space for the basis functions 
N = zeros(d+1,1);                               %   double *N = (double*) mxMalloc((d+1)*sizeof(double)); 
                                                % 
                                                %   // for each parametric point i 
for col=1:nu                                    %   for (col = 0; col < nu; col++) { 
                                                %     // find the span of u[col] 
    s = findspan(nc-1, d, u(col), k);           %     s = findspan(nc-1, d, u[col], k); 
    N = basisfun(s,u(col),d,k);                 %     basisfun(s, u[col], d, k, N); 
                                                % 
    tmp1 = s - d + 1;                           %     tmp1 = s - d; 
    for row=1:mc                                %     for (row = 0; row < mc; row++)  { 
        tmp2 = 0;                               %       tmp2 = 0.0; 
        for i=0:d                               %       for (i = 0; i <= d; i++) 
           tmp2 = tmp2 + N(i+1)*c(row,tmp1+i);  % 	tmp2 += N[i] * ctrl[tmp1+i][row]; 
        end                                     % 
        p(row,col) = tmp2;                      %       pnt[col][row] = tmp2; 
    end                                         %     } 
end                                             %   } 
                                                % 
                                                %   mxFree(N); 
                                                %   freevec2mat(pnt); 
                                                %   freevec2mat(ctrl); 
                                                % 
                                                %   return ierr; 
                                                %   } 
end


function s = findspan(n,p,u,U)                 
% FINDSPAN  Find the span of a B-Spline knot vector at a parametric point 
% ------------------------------------------------------------------------- 
% ADAPTATION of FINDSPAN from C 
% ------------------------------------------------------------------------- 
% 
% Calling Sequence: 
%  
%   s = findspan(n,p,u,U) 
%  
%  INPUT: 
%  
%    n - number of control points - 1 
%    p - spline degree 
%    u - parametric point 
%    U - knot sequence 
%  
%  RETURN: 
%  
%    s - knot span 
%  
%  Algorithm A2.1 from 'The NURBS BOOK' pg68 
                                                 
                                                % int findspan(int n, int p, double u, double *U) { 
                                                 
                                                %   int low, high, mid;                                                 
                                                %   // special case 
if (u==U(n+2)), s=n; return,  end               %   if (u == U[n+1]) return(n); 
                                                % 
                                                %   // do binary search 
low = p;                                        %   low = p; 
high = n + 1;                                   %   high = n + 1; 
mid = floor((low + high) / 2);                  %   mid = (low + high) / 2; 
while (u < U(mid+1) || u >= U(mid+2))           %   while (u < U[mid] || u >= U[mid+1])  { 
    if (u < U(mid+1))                           %     if (u < U[mid]) 
        high = mid;                             %       high = mid; 
    else                                        %     else 
        low = mid;                              %       low = mid;                   
    end  
    mid = floor((low + high) / 2);              %     mid = (low + high) / 2; 
end                                             %   } 
                                                % 
s = mid;                                        %   return(mid); 
                                                %   } 
end


function N = basisfun(i,u,p,U)                 
% BASISFUN  Basis function for B-Spline 
% ------------------------------------------------------------------------- 
% ADAPTATION of BASISFUN from C Routine 
% ------------------------------------------------------------------------- 
% 
% Calling Sequence: 
%  
%   N = basisfun(i,u,p,U) 
%    
%    INPUT: 
%    
%      i - knot span  ( from FindSpan() ) 
%      u - parametric point 
%      p - spline degree 
%      U - knot sequence 
%    
%    OUTPUT: 
%    
%      N - Basis functions vector[p+1] 
%    
%    Algorithm A2.2 from 'The NURBS BOOK' pg70. 
                                                 
                                                  %   void basisfun(int i, double u, int p, double *U, double *N) { 
                                                  %   int j,r; 
                                                  %   double saved, temp; 
i = i + 1; 
                                                  %   // work space 
left = zeros(p+1,1);                              %   double *left  = (double*) mxMalloc((p+1)*sizeof(double)); 
right = zeros(p+1,1);                             %   double *right = (double*) mxMalloc((p+1)*sizeof(double)); 
                                                
N(1) = 1;                                         %   N[0] = 1.0; 
for j=1:p                                         %   for (j = 1; j <= p; j++) { 
    left(j+1) = u - U(i+1-j);                     %   left[j]  = u - U[i+1-j]; 
    right(j+1) = U(i+j) - u;                      %   right[j] = U[i+j] - u; 
    saved = 0;                                    %   saved = 0.0; 
 
    for r=0:j-1                                   %   for (r = 0; r < j; r++) { 
        temp = N(r+1)/(right(r+2) + left(j-r+1)); %   temp = N[r] / (right[r+1] + left[j-r]); 
        N(r+1) = saved + right(r+2)*temp;         %   N[r] = saved + right[r+1] * temp; 
        saved = left(j-r+1)*temp;                 %   saved = left[j-r] * temp; 
    end                                           %   } 
 
    N(j+1) = saved;                               %   N[j] = saved; 
end                                               %   } 
   
                                                  %   mxFree(left); 
                                                  %   mxFree(right); 
                                                  %   } 
end