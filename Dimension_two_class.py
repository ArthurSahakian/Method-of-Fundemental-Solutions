# -*- coding: utf-8 -*-
"""
Created on Sun Apr  6 20:47:15 2025

@author: sahakian.a
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import scipy


class dimension_2():
    
    '''
    This class computes the eigenvalues and eigenvectors of the biharmonic Steklov
    
    eigenvalue problem on different planar shapes, using the method of fundamental 
    
    solutions.

    Parameters
    ----------
    N : float or integer
    Number of points considered on the boundary of the domain Omega. 
    beta : float or integer
    The distance between the boundary of Omega, and the artificial boundary of Omega.
    Returns
    -------
    An instance of the class.

    '''
    
    def __init__(self,N,beta=0.4):
        self.N=N
        self.h=2*np.pi/N
        self.beta=beta
    
    def fourrier_fit(self,params):
        '''
        
        Computes the support function p, as well as its derivative, based on the fourrier
        
        aproximations, where the "params" variable is the the fourrier coefficient of the
        
        support function.
        

        Parameters
        ----------
        params : numpy array or list
            The parameters of the Fourrier approximation of the support function
        

        '''
        n = int(np.size(params)/2) + 1
        main_array_cos = np.arange(n)
        main_array_sin = np.arange(1, n)
        
        def n_func_p(x):
            result_cos = np.cos(x[:, np.newaxis] * main_array_cos)
            result_sin = np.sin(x[:, np.newaxis] *main_array_sin)
            res= np.sum( result_cos* params[0:n],axis=1) + np.sum( result_sin*params[n:],axis=1)
            return res
        
        def n_p_prime(x):
            result_cos = -np.sin(x[:, np.newaxis] * main_array_cos)
            result_sin = np.cos(x[:, np.newaxis] *main_array_sin)
            res= np.sum( result_cos* params[0:n] *main_array_cos ,axis=1) + np.sum( result_sin*params[n:] * main_array_sin,axis=1)
            return res
        
        def n_sum_p_p_second(x):
           result_cos = np.cos(x[:, np.newaxis] * main_array_cos)
           result_sin = np.sin(x[:, np.newaxis] *main_array_sin)
           res= np.sum( result_cos* params[0:n] * (1- main_array_cos**2) ,axis=1) + np.sum( result_sin*params[n:] * (1-main_array_sin**2),axis=1)
           return np.abs(res)
       
        self.p_prime= n_p_prime
        self.sum_p_p_second = n_sum_p_p_second
        self.func_p= n_func_p
        self.params=params
        
    def angles_gen(self):
        '''
        Generates equidistant angles according to the number of points considered.
        
        These angles will be useful to represent the points in polar coordinates.

        Returns
        -------
        angles : float
            N equidistant angles between 0 and 2 PI.

        '''
        try:
            params = self.params
        except AttributeError:
            raise Exception("You need to call the method fourrier_fit with appropriate Fourrier parameters that describe the support function of the surface of your domain.")

        n = int(np.size(params)/2) + 1
        main_array_cos = np.arange(1,n)
        main_array_sin = np.arange(1, n)
        
        def n_func_arc(angle):
            first=angle*params[0]
            result_cos = np.sin(angle * main_array_cos)
            result_sin = -np.cos(angle *main_array_sin)
            res= np.sum( result_cos* params[1:n] *(1-main_array_cos**2)/main_array_cos ) + np.sum( result_sin*params[n:] *(1-main_array_sin**2)/main_array_sin  )
            return res + first

        ### here angles with arc length
        try:
            object_total_arc= n_func_arc(2*np.pi) - n_func_arc(0)
            angles=np.zeros([self.N])
            for i in range(0,self.N-1):
                local_func = lambda t : n_func_arc(t) - n_func_arc(angles[i]) - object_total_arc/(self.N) 
                res = scipy.optimize.root_scalar(local_func,bracket=[angles[i],2*np.pi],xtol=1e-12)
                angles[i+1] = res.root
        except:
            print("Specialized angle choser failed. Angles were chosen equidistantly inside [0,2 pi]. This will affect performance.")
            angles= np.array([i*2*np.pi/self.N for i in range(self.N)])
        self.angles=angles
        return angles
        
    def X_given_angle(self,t):
        """
        This function takes an array of angles between 0 and 2 pi as input and returns 
        the associated point on the boundary.

        Parameters
        ----------
        t : numpy array
            A 1D numpy array containing a list of angles.

        Returns
        -------
        X : numpy array
            A 2D numpy array containing the associated points to the selected angles.

        """
        P= self.func_p(t)
        Pprime= self.p_prime(t)
        Cos= np.cos(t)
        Sin= np.sin(t)
        X=np.zeros([len(t),2])
        X[:,0]= Cos*P-Pprime*Sin
        X[:,1]=Sin*P + Cos* Pprime
        return X
    
    def normal_given_angle(self,t):
        """
        This function takes an array of angles between 0 and 2 pi as input and returns 
        the normals of the associated point on the boundary.
        
        Parameters
        ----------
        t : numpy array
            A 1D numpy array containing a list of angles.

        Returns
        -------
        normals : numpy array
            A 2D numpy array containing the associated normals to the selected angles.

        """
        normals = np.zeros([len(t),2])
        normals[:,0]= np.cos(t)
        normals[:,1]= np.sin(t)
        
        return normals
    
    def phi(self,x):
        '''
        The fundamental solution function of the Bilaplacien.

        Parameters
        ----------
        x : float
            a norm

        Returns
        -------
        TYPE float
            The value of the function phi at the the point x.

        '''
        return (x**2)*np.log(x)
    
    def matrix_calcul_steklov(self,X,normals):
        """
        Computes the matrices Z and W of the generalized eigenvalue problem Zv=d Wv, associated with the Biharmonic Steklov problem.

        Parameters
        ----------
        X : 2D numpy array
            The list of points describing the boundary of the domain.
        normals : 2D numpy array
            The normalized normals at the points X.

        Returns
        -------
        Z : 2D numpy array
            Matrix on the left hand side of the equation.
        W : 2D numpy array
            Matrix on the right hand side of the equation.

        """
        
        normalized_normals=normals
        normalized_normals_Y= normals
        Y= X + self.beta*normalized_normals_Y
        Z=np.zeros([2*self.N,2*self.N],dtype=np.float64)
        W=np.zeros([2*self.N,2*self.N],dtype=np.float64)
        di = X[:, np.newaxis, :] - Y[np.newaxis, :, :]
        norm_di= np.linalg.norm(di,axis=2)
        di_scalar_n = np.sum(normalized_normals[:, np.newaxis, :] * di, axis=-1)
        di_scalar_ni = np.sum( di* normalized_normals, axis=2)
        n_scalar_ni = np.sum(normalized_normals[:,np.newaxis,:]*normalized_normals[np.newaxis,:,:],axis=2)
        Z[:self.N,:self.N]= self.phi(norm_di)
        Z[:self.N,self.N:]= (2 * np.log(norm_di) +1) * di_scalar_ni
        Z[self.N:,:self.N]= 4*( np.log(norm_di) +1 )
        Z[self.N:,self.N:] = 4 * di_scalar_ni / (norm_di)**2
        W[self.N:,:self.N]= -di_scalar_n * (2 * np.log(norm_di) +1) 
        W[self.N:,self.N:]= - 2*di_scalar_n*di_scalar_ni/ (norm_di**2) - (2 * np.log(norm_di) +1) * n_scalar_ni
        self.Z=Z
        self.W=W
        self.Y=Y
        self.A= Z[:self.N,:self.N]
        return Z, W
    
    def matrix_calcul_variant(self,X,normals):
        """
        Computes the matrices Z and W of the generalized eigenvalue problem Zv=d Wv, associated with the Biharmonic Steklov variant problem.

        Parameters
        ----------
        X : 2D numpy array
            The list of points describing the boundary of the domain.
        normals : 2D numpy array
            The normalized normals at the points X.

        Returns
        -------
        Z : 2D numpy array
            Matrix on the left hand side of the equation.
        W : 2D numpy array
            Matrix on the right hand side of the equation.

        """
        
        normalized_normals=normals
        normalized_normals_Y= normals
        Y= X + self.beta*normalized_normals_Y
        Z=np.zeros([2*self.N,2*self.N],dtype=np.float64)
        W=np.zeros([2*self.N,2*self.N],dtype=np.float64)
        Z=np.zeros([2*self.N,2*self.N])
        W=np.zeros([2*self.N,2*self.N])
        di = X[:, np.newaxis, :] - Y[np.newaxis, :, :]
        norm_di= np.linalg.norm(di,axis=2)
        di_scalar_n = np.sum(normalized_normals[:, np.newaxis, :] * di, axis=-1)
        di_scalar_ni = np.sum( di* normalized_normals, axis=2)
        n_scalar_ni = np.sum(normalized_normals[:,np.newaxis,:]*normalized_normals[np.newaxis,:,:],axis=2)
        Z[:self.N,:self.N]= self.phi(norm_di)
        Z[:self.N,self.N:]= (2 * np.log(norm_di) +1) * di_scalar_ni
        Z[self.N:,:self.N]= 4*( np.log(norm_di) +1 )
        Z[self.N:,self.N:] = 4 * di_scalar_ni / (norm_di)**2 
        W[self.N:,:self.N]= -di_scalar_n * (2 * np.log(norm_di) +1) 
        W[self.N:,self.N:]= - 2*di_scalar_n*di_scalar_ni/ (norm_di**2) - (2 * np.log(norm_di) +1) * n_scalar_ni
        angles = self.angles_gen()
        Pp = 1/ self.sum_p_p_second(angles)
        modify_pp= Pp[np.newaxis,:].T
        Z[self.N:,:self.N] +=  W[self.N:,:self.N] * modify_pp
        Z[self.N:,self.N:] += W[self.N:,self.N:] * modify_pp
        self.Z=Z
        self.W=W
        self.A= Z[:self.N,:self.N]
        
        return Z,W 
    
        
    def plot_shape(self, scatter=False, ann=False, include_art_boundary=False):
        """
        Plots the shape described by the domain. Optionally includes scatter plot, annotations, and artificial boundary.
    
        Parameters
        ----------
        scatter : bool, optional
            If True, uses scatter plot instead of a connected line plot. Default is False.
        ann : bool, optional
            If True, annotates the points with their indices. Default is False.
        include_art_boundary : bool, optional
            If True, plots the artificial boundary in red. Default is False.
        """
    
        # Ensure self.X and normals are defined
        if not hasattr(self, 'X') or self.X is None:
            angles = self.angles_gen()
            self.X = self.X_given_angle(angles)
            normals = self.normal_given_angle(angles)
        else:
            angles = self.angles_gen()
            normals = self.normal_given_angle(angles)
    
        Y = self.X + self.beta * normals
    
        # Prepare closed-loop shape
        Z = np.vstack([self.X, self.X[0]])
    
        fig, ax = plt.subplots()
    
        # Plot original shape
        if scatter:
            ax.scatter(Z[:, 0], Z[:, 1], color='b', label='Original Shape')
        else:
            ax.plot(Z[:, 0], Z[:, 1], color='b', label='Original Shape')
    
        # Annotations
        if ann:
            for i, (x, y) in enumerate(Z):
                ax.annotate(str(i), (x, y))
    
        # Plot artificial boundary
        if include_art_boundary:
            W = np.vstack([Y, Y[0]])
            if scatter:
                ax.scatter(W[:, 0], W[:, 1], color='r', label='Artificial Boundary')
            else:
                ax.plot(W[:, 0], W[:, 1], color='r', label='Artificial Boundary')
    
        ax.set_aspect('equal')
        ax.legend()
        plt.show()


            
    def real_eigenvalue(self,i=0,method="direct",problem="steklov"):
        '''
        Compute the approximate eigenvalues by the method of the fundamental solution 
        S means we have a square generalized eigenvalue problem.

        Parameters
        ----------
        i : integer, optional
           The rank of the eigenvalue. The default is 0.
        method : string, optional
            Method of computation of the eigenvalue problem. Options are "direct",
            "inverse" and "lu". The default is direct.
        problem : string, optional
            Controls the problem solved by this function. The options are "steklov" (original problem) and "variant" (the variant problem).
        Returns
        -------
        first_eigen : float
            The eigenvalue of rank i.
        positive : list of floats
            All the eigenvalues of the generalized problem.
        '''
        
        ##Generate angles, points and normals
        angles=self.angles_gen()
        X=self.X_given_angle(angles)
        normals= self.normal_given_angle(angles)
        
        ##Compute the matrices according to problem
        if problem == "steklov":
            Z,W =  self.matrix_calcul_steklov(X,normals)
            
        elif problem == "variant":
            Z,W =  self.matrix_calcul_variant(X,normals)
        
        if method == "direct":
            eig=scipy.linalg.eig(self.Z,-self.W)
            
        elif method == "inverse":
            inverse= np.linalg.inv( np.copy(self.Z))
            multi = np.matmul( -self.W,inverse)
            eig = np.linalg.eig(multi)
        
        elif method == "lu":
            l, u = scipy.linalg.lu(self.Z.copy(),permute_l=True)
            l_1 , u_1 = scipy.linalg.lu(-self.W.copy(),permute_l=True)
            truc_muche= np.linalg.inv(l) @ l_1 @ u_1 @  np.linalg.inv(u)
            eig= scipy.linalg.eig(truc_muche)
            
        else:
            Q, R = np.linalg.qr(-self.W)
            eig=scipy.linalg.eig(np.matmul(Q.T,self.Z) ,R)
        eigens=eig[0]
        real_eigens=eigens[np.isreal(eigens)]
        if method =="inverse":
            positive= np.real(np.sort(1/real_eigens[real_eigens>0]))
        elif method =="lu":
            positive= np.real(np.sort(1/real_eigens[real_eigens>0]))
        else:
            positive=np.real(np.sort(real_eigens[real_eigens>0]))
        first_eigen=positive[i]
        smallest_eigenvalue_index =  eigens== first_eigen
        eigenvector= np.real(eig[1] [:,smallest_eigenvalue_index])
        self.eigenvalue=first_eigen
        self.eigenvector=eigenvector    
        return first_eigen, positive
    
        
    