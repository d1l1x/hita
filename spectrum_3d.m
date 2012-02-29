%% Clear complete workspace
% <latex>
% If you start a new Matlab project it is always a goot idea
% to clear the complete workspace and command window.
% Since Matlab does not automatically open new figure windows
% each time a plot is invoked we should close
% all remaining figure windows. Also very important is to tell
% Matlab where to find the functions we will use in our project.
% This is done by adding the functions folder with
% \lstinline!path('./functions',path)! to the Matlab search
% path.
% </latex>
%
display('Clear workspace ...')
path('./functions',path) % add directory the Matlab path
close all % close all figures
clear all % clear workspace
clc % clear command window

datadir='data/stuttgart/'; % data directory
flag='3D';
%set(0,'DefaultFigureWindowStyle','docked')
%% Read data files
%
% <latex>
% During the evaluation of the function \verb|ReadData| all velocity files
% necessary for the calculation of the spectrum and the correlation
% coefficients are read. In addition the
% import operations are enclosed in a \lstinline!tic;!$\ldots$
% \lstinline!;toc! block which
% measures the time needed for reading the ASCII data.
% Although the ASCII
% data format ist not the best choice in terms of speed and size, we
% will use it since other
% methodologies require additional knowledge of data processing. Just for
% your information a very famous and highly portable data format is  
% \href{http://www.hdfgroup.org/HDF5/}{hdf5}. It is a software library that
% is available for a range of computer platforms, from laptops to massively
% parallel systems and implements a high-level API (Application 
% programming interface) with C, \C++, Fortran 90,
% and Java interfaces. Besides its hierarchical
% structure it is highly optimized for parallel I/O operations and can be
% read by nearly all data processing/visualization tools.
% </latex>
%
display('Read data ...')
[uvel,vvel,wvel,time_read,dim] = ReadData(datadir,flag,...
                                                 'uvel',...
                                                 'vvel',...
                                                 'wvel');
%% 
% <latex>
% \lstinputlisting{../functions/ReadData.m}
% </latex>
%% Set necessary parameters
% 
% <latex>
% For further computations we have to to define some parmeters of the
% DNS simulation.
% </latex>
%
display('Set parameters ...')
[u,v,w,Lx,nu]=Params(uvel,vvel,wvel,dim);
%% 
% <latex>
% \lstinputlisting{../functions/Params.m}
% </latex>
%% Compute 3D spectrum
% 
% <latex>
% The core of the code is contained in the function
% \lstinline!PowerSpec!. It computes the three dimensional energy spectrum
% from the given velocity fields, obtained from a direct numerical
% simulation. Although the theoretical analysis is 
% relatively demanding compared to one dimensional spectra its worth
% investing the effort.
% The theory of one dimensional spectra relies
% on the assumption that the propagation of spectral waves ($\kappa_1$)
% is in the direction of the observed velocity fields or to say it differently one
% dimensional spectra and correlation functions are Fourier transform pairs.
% The theory of correlation functions will be discussed in section \ref{sec:correlation}.
% A key drawback of this theory is that the calculated spectrum has
% contributions from all wavenumbers $\boldsymbol\kappa$, so that the
% magnitude of $\boldsymbol\kappa$ can be appreciably larger than 
% $\kappa_1$. This phenomenon is called aliasing.
% In order to avoid these aliasing effects is also possible to produce correlations that
% involve all possible directions. The three dimensional Fourier
% transformation of such a correlation produces a spectrum that not only
% depends on a single wavenumber but on the wavenumber vector $\kappa_i$.
% Though the directional information contained in $\kappa_i$ eliminates the
% aliasing problem the complexity makes a physical reasoning impossible.
% For homogeneous isotropic turbulence the situation can be considerably
% simplified. From the knowledge that the velocity field is isotropic it can be
% shown that the velocity spectrum tensor is fully determined by
%   \begin{equation}
%       \label{eq:iso_tensor}
%       \Phi_{ij}(\boldsymbol\kappa) = A(\kappa)\delta_{ij}+B(\kappa)\kappa_i\kappa_j,
%   \end{equation}
% where $A(\kappa)$ and $B(\kappa)$ are arbitrary scalar functions. Since we assume
% incompressible fluids (mathematically expressed by $\nabla\cdot u=0$ or $\kappa_iu_i=0$
% the following condition holds
%   \begin{equation}
%       \kappa_i\,\Phi_{ij}(\boldsymbol\kappa)=0.
%   \end{equation}
% It can be shown that this yields a relation between $A$ and $B$ by means of
%   \begin{equation}
%       \label{eq:rel_AB}
%   	B(|\boldsymbol\kappa|)=-\frac{A(|\boldsymbol\kappa|)}{(|\boldsymbol\kappa|)^2}	
%   \end{equation}
% In the end this gives a relation between the three dimensional energy spectrum
% function $E(|\boldsymbol\kappa|)$ and the velocity spectrum tensor $\Phi_{ij}$.
%   \begin{equation}
%   	\Phi_{ij}=\frac{E(|\boldsymbol\kappa|)}{4\pi\,(|\boldsymbol\kappa|)^2}\left(\delta_{ij}
%   	-\frac{\kappa_i\kappa_j}{(|\boldsymbol\kappa|)^2}\right)	
%   \end{equation}
% The question is now how the remaining variable ($A$ or $B$) can be determined. Regarding the turbulent kinetic
% energy we know that
%   \begin{equation}
%       \label{eq:exp_for_k}
%       k=\int\limits_{-\infty}^{\infty}E(|\boldsymbol\kappa|)\,\mathrm{d}k
%       =\sum\limits_{\boldsymbol\kappa}E(\boldsymbol\kappa)
%       =\sum\limits_{\boldsymbol\kappa}\frac{1}{2}\left<u^{*}(\boldsymbol\kappa)\,u(\boldsymbol\kappa)\right>
%       =\iiint\limits_{-\infty}^{\infty}\frac{1}{2}\Phi_{ii}(\boldsymbol\kappa)\,\mathrm{d}\boldsymbol\kappa.
%   \end{equation}
% Comparing the second and last expression we get
%   \begin{equation}
%       E(|\boldsymbol\kappa|)=\oiint\frac{1}{2}\,\Phi_{ii}(\boldsymbol\kappa)\,\mathrm{d}S(\kappa).
%   \end{equation}
%   \begin{figure}[t!]
%       \centering
%       \includegraphics[scale=1]{shell_integration}
%       \label{fig:shell_int}
%       \caption{Illustration of the two dimensional shell integration}
%   \end{figure}
% This integral can be solved analytically by utilizing again the assumption of isotropy.
% For these kind of flows the energy spectrum function can be regarded as the sum of kinetic energy
% (in wave number space) on different energy levels. Each of these energy levels is denoted by a spherical
% shell in wave number space. Since the surface of a sphere is completely determined by its radius the
% surface integral can be solved analytically. The idea of this integration is illustrated
% in Fig. \ref{fig:shell_int}.
% As a result of this one gets
% \begin{figure}[t!]
%       \includegraphics[scale=0.7]{spectrum}
%       \caption{Computed spectrum}
%       \label{fig:shell_int}
% \end{figure}
%   \begin{equation}
%       E(|\boldsymbol\kappa|)=\oiint\frac{1}{2}\,\Phi_{ii}(\boldsymbol\kappa)\,\mathrm{d}S(\kappa)
%       =4\pi(|\boldsymbol\kappa|)^2\,\Phi_{ii}(|\boldsymbol\kappa|).
%   \end{equation}
% Introducing this relation to equations \eqref{eq:iso_tensor} and
% \eqref{eq:rel_AB} one arrives at an expression for the variable $B$.
%   \begin{equation}
%       B=-\frac{E(|\boldsymbol\kappa|)}{4\pi(|\boldsymbol\kappa|)^2} 
%   \end{equation}
% Together with the approximation of the integral of $\Phi$
% (equation \eqref{eq:exp_for_k})
%   \begin{equation}
%       \label{eq:Phi_approx}
%       \iiint\limits_{-\infty}^{\infty}\frac{1}{2}\Phi_{ii}(\boldsymbol\kappa)\,\mathrm{d}\boldsymbol\kappa
%       \approx\frac{1}{2}\sum\limits_{\boldsymbol\kappa}\Phi_{ii}(\boldsymbol\kappa)
%       \,\Delta\kappa_x\,\Delta\kappa_y\,\Delta\kappa_z,
%   \end{equation}
% where $\Delta\kappa$ refers to the step size in wave number space,
% the final expression of the three dimensional discrete energy spectrum can be derived.
%   \begin{equation}
%       E(|\boldsymbol\kappa|)=2\pi(|\boldsymbol\kappa|)^2\frac{\left<u^{*}
%       (\boldsymbol\kappa)\,u(\boldsymbol\kappa)\right>}{(\Delta\kappa)^3}
%   \end{equation}
% The calling sequence for the computation of the energy spectrum reads
% </latex>
display('Compute spectrum...')
[spectrum,k,bin_counter,time_spec,R,dx] = PowerSpec(u,...
                                                    v,...
                                                    w,...
                                                  Lx,dim);
%% 
% <latex>
% \lstinputlisting{../functions/PowerSpec.m}
% </latex>
%% Compute dissipation and turbulent kinetic energy
%
% <latex>
%   The function \lstinline!SpecProp! calculates the kinetic energy both
%   from the velocities and the previously computed spectrum. The latter
%   one is calculated by using relation \eqref{eq:exp_for_k} and
%   incorporating equation \eqref{eq:Phi_approx}
%   \begin{equation}
%       k\approx\sum\limits_{\boldsymbol\kappa}\frac{E(|\boldsymbol\kappa|)}{4\pi\,
%       (|\boldsymbol\kappa|)^2}\,\Delta\kappa_x\,\Delta\kappa_y\,
%       \Delta\kappa_z
%   \end{equation}
% Knowing that the underlying flow field is isotropic this triple sum may be
% reduced to a single one by
%   \begin{equation}
%       k\approx\sum\limits_{\boldsymbol\kappa}\frac{E(|\boldsymbol\kappa|)}{4\pi\,
%       (|\boldsymbol\kappa|)^2}\,\Delta\kappa_x\,\Delta\kappa_y\,
%       \Delta\kappa_z=
%       \sum\limits_{|\boldsymbol\kappa|}E(|\boldsymbol\kappa|)\cdot\Delta|\boldsymbol\kappa|
%   \end{equation}
%   A second integral, approximated in the same manner, gives the value
%   of the Dissipation
%   \begin{align}
%       \label{eq:Phi_approx}
%       \varepsilon=\iiint\limits_{-\infty}^{\infty}2(|\boldsymbol\kappa|)^2\nu\,\frac{1}{2}\Phi_{ii}(\boldsymbol\kappa)\,\mathrm{d}\boldsymbol\kappa
%       &\displaystyle\approx\frac{1}{2}\sum\limits_{\boldsymbol\kappa}2(|\boldsymbol\kappa|)^2\nu\,\Phi_{ii}(\boldsymbol\kappa)
%       \,\Delta\kappa_x\,\Delta\kappa_y\,\Delta\kappa_z\\
%       &\displaystyle\approx\sum\limits_{|\boldsymbol\kappa|}2(|\boldsymbol\kappa|)^2
%       \nu\,E(|\boldsymbol\kappa|)\,\Delta(|\boldsymbol\kappa|)
%   \end{align}
%   where $\nu$ refers to the kinematic viscosity.
%   The calling sequence reads
% </latex>
display('Compute kinetic energy...')
[Dissipation,kin_Sp,kin_Ph,kin_E,up] = SpecProp(spectrum,k,...
                                                nu,u,v,w,...
                                                dim,dx);
%% 
% <latex>
% \lstinputlisting{../functions/SpecProp.m}
% </latex>
%% Kolmogorov properties
%
% <latex>
% According to the Kolomogorov hypotheses the length scale $\eta$,
% characteristic velocity $u_{\eta}$ and the characteristic time scale
% $\tau$ of the smallest swirls in the flow are computed
% within the function \lstinline!KolmoScale!. From a dimensionality
% analysis Kolmogorov derived
%   \begin{eqnarray}
%       \eta&=\displaystyle\left(\frac{\nu^3}{\varepsilon}\right)^{1/4},\\
%       u_{\eta}&=\left(\varepsilon\,\nu\right)^{1/4},\\
%       \tau_{\eta}&=\displaystyle\left(\frac{\nu}{\varepsilon}\right)^{1/2}.
%   \end{eqnarray}
% For further reading concerning his theory it is refered to
% \citet{Pope:2000tp}, \citet{Hinze:1975tb} and \citet{Tennekes:1972vb}.
% </latex>
display('Compute Kolmogorov scales...')
[eta,u_eta,tau]=KolmoScale(nu,Dissipation);
%% 
% <latex>
% \lstinputlisting{../functions/KolmoScale.m}
% </latex>
%% Compute correlations
%
% <latex>
% \label{sec:correlation}
% From a general perspective correlation functions are a measure of how
% much two physical quantities are connected. So how is this helpful for
% the analysis of turbulent flows? For seemingly chaotic and random
% processes it would be beneficial if we had a measure of how the velocity
% at point $A$ is influenced by the velocity at point $B$. A maybe more
% intuitive quantity that can be calculated from the correlation functions
% is the integral length scale which gives a measure of the largest
% eddies in the flow. In fluid dynamics one generally differentiates
% between two forms of correlation functions, the longitudinal and the
% transversal or lateral correlation function. The difference between both forms is
% illustrated in figure \ref{fig:correlations}.
% \begin{figure}[t!]
%         \centering
%         \subfloat[Longitudinal correlation][Longitudinal correlation]{
%         \includegraphics[scale=1]{long_corr.eps}
%         \label{fig:long_corr}
%         }
%         \hfill
%         \subfloat[Transversal correlation][Transversal correlation]{
%         \includegraphics[scale=1]{trans_corr.eps}
%         \label{fig:trans_corr}
%         }
% \caption{Illustration of different correlation functions}
% \label{fig:correlations}
% \end{figure}
% In general the correlation between two components of an isotropic
% homogeneous velocity field is expressed by
%   \begin{equation}
%       R_{ij} = \left<u_i(\mathbf{x}+r)\,u_j(\mathbf{x})\right>.
%   \end{equation}
% The correlation coefficients are computed by normalizing
% $R_{ij}$ with the square root of the product of the two variances
% $\sigma_i^2$ and $\sigma_j^2$.
%   \begin{equation}
%       r_{ij} = \frac{\left<u_i(\mathbf{x}+\mathbf{r})\,u_j(\mathbf{x})\right>}
%                   {\sqrt{\sigma_i^2\,\sigma_j^2}}
%   \end{equation}
% As illustrated in figure \ref{fig:correlations} the longitudinal and lateral
% correlation depend on the direction of $\mathbf{r}$, i.e.
%   \begin{eqnarray}
%       f(r) &= \displaystyle\frac{\left<u_1(\mathbf{x}+r\mathbf{e}_1)\,u_1(\mathbf{x})\right>}
%                   {\sigma_1}\\
%       g(r) &= \displaystyle\frac{\left<u_2(\mathbf{x}+r\mathbf{e}_1)\,u_2(\mathbf{x})\right>}
%                   {\sigma_2}
%   \end{eqnarray}
% From several theoretical analysis it is well known that correlations can
% be efficiently computed by means of multiplying the Fourier transform of
% a quantity with its complex conjugate. Using the FFT approach this gives
% an enormeous speed advantage.
%   \begin{equation}
%       R_{ij} =
%       \mathfrak{F}^{-1}\left\{\mathfrak{F}\left\{u_i\right\}^*\cdot\mathfrak{F}\left\{u_j\right\}\right\}
%   \end{equation}
% </latex>
display('Compute Correlations...')
[R11,R22,r,R1,R2,R3]=Correlation(u,v,w,Lx,dim);
%% 
% <latex>
% The content of \verb|Correlation| reads
% \lstinputlisting{../functions/Correlation.m}
% </latex>
%%
% 
% <latex>
% \bibliographystyle{NTFD-bibstyle}
% \bibliography{bibliography}
% </latex>