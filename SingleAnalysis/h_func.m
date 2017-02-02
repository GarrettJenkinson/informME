%   informME: An information-theoretic pipeline for WGBS data
%   Copyright (C) 2017, Garrett Jenkinson (jenkinson@jhu.edu)
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software Foundation,
%   Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
%   or see <http://www.gnu.org/licenses/>.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  informME: Information-Theoretic Analysis of Methylation  %%%%%%%%
%%%%%%%%                        h_func.m                           %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
%%%%%%%%               Last Modified: 11/30/2016                   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function computes the (log2-based) entropies of a collection 
% X(q,r), q = 1,2,...,Q, r = 1,2,...,R of QxR binary random variables 
% with corresponding probabilities p(q,r) and 1-p(q,r). 
%
% USAGE: ENTR = h_func(P)
%
% INPUT:
%
% P
%         A QxR matrix of the probabilities p(q,r), q = 1,2,...,Q, 
%         r = 1,2,...,R. 
%
% OUTPUT:
%
% ENTR    
%         A QxR matrix with its (q,r) element being the entropy of 
%         the binary random variable X(q,r).
%       

function ENTR = h_func(P)

ENTR          = zeros(size(P));
ENTR(P>0&P<1) = -P(P>0&P<1).*log2(P(P>0&P<1))...
                -(1-P(P>0&P<1)).*log2(1-P(P>0&P<1));