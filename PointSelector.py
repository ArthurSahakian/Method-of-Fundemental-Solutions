# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 10:49:03 2024

@author: sahakian.a
"""

import numpy as np
import scipy
import matplotlib.pyplot as plt

class point_selector_dimension_two:
    
    def __init__(self):
        pass 
    
    def fit_fourier(self,N,fourier_coeff,domain):
        """
        Fit a Fourier series to the given data.
    
        Parameters:
        - N: Number of terms in the Fourier series.
        - fourier_coeff: Fourier coefficients for the series.
        - domain: Domain of the data.
    
        Usage example:
        obj.fit_fourier(N=10, fourier_coeff=[...], domain=[start_value, end_value])
    
        This method sets the parameters of a Fourier fit in the class, including the number of terms (N),
        the Fourier coefficients (fourier_coeff), and the domain of the data (domain).
    
        Note: Make sure to provide valid values for N, fourier_coeff, and domain.
    
        After calling this method, you may have additional methods or attributes to perform operations or evaluations based on the Fourier fit.

       """
        self.N=N
        self.fourier_coeff = fourier_coeff
        self.domain = domain
        
    def fourier_func_initialization(self):
        """
        Initialize Fourier function based on the provided coefficients.
    
        This function initializes a Fourier function based on the given Fourier coefficients.
        It calculates the number of coefficients, defines arrays for cosine and sine terms,
        and creates functions to evaluate the Fourier series and its derivatives.
    
        Parameters:
        - None (Uses attributes of the class)
    
        Returns:
        - None (Assigns functions to class attributes)
       """
        n = int(np.size(self.fourier_coeff)/2) + 1
        main_array_cos = np.arange(n)
        main_array_sin = np.arange(1, n)
        def func_p(x):
            result_cos = np.cos(x[:, np.newaxis] * main_array_cos)
            result_sin = np.sin(x[:, np.newaxis] *main_array_sin)
            res= np.sum( result_cos* self.fourier_coeff[0:n],axis=1) + np.sum( result_sin*self.fourier_coeff[n:],axis=1)
            return res
        def func_p_prime(x):
            result_cos = -np.sin(x[:, np.newaxis] * main_array_cos)
            result_sin = np.cos(x[:, np.newaxis] *main_array_sin)
            res= np.sum( result_cos* self.fourier_coeff[0:n] *main_array_cos ,axis=1) + np.sum( result_sin*self.fourier_coeff[n:] * main_array_sin,axis=1)
            return res
        def func_p_second(x):
           result_cos = np.cos(x[:, np.newaxis] * main_array_cos)
           result_sin = np.sin(x[:, np.newaxis] *main_array_sin)
           res= np.sum( result_cos* self.fourier_coeff[0:n] * (1- main_array_cos**2) ,axis=1) + np.sum( result_sin*self.fourier_coeff[n:] * (1-main_array_sin**2),axis=1)
           return res
        def func_arc(x):
            first_comp=x*self.fourier_coeff[0]
            main_array_cos = np.arange(1,n)
            result_cos = np.sin(x * main_array_cos)
            result_sin = -np.cos(x *main_array_sin)
            second_comp= np.sum( result_cos* self.fourier_coeff[1:n] *(1-main_array_cos**2)/main_array_cos ) + np.sum( result_sin*self.fourier_coeff[n:] *(1-main_array_sin**2)/main_array_sin  )
            return first_comp + second_comp
        self.func_p= func_p
        self.func_p_prime= func_p_prime
        self.func_p_second= func_p_second
        self.func_arc=func_arc
        
    def angles_gen_arc(self,p,domain=(0,2*np.pi)):
        """
        Generates equidistant angles on an arc, given the parametric equation of its boundary.
        By default, the domain of definition of the function is [0,2pi].

        Parameters
        ----------
        p : array
            fourrier coefficients of the boundary function.
        domain : tuple or list
                 domain of definition of the function.

        Returns
        -------
        angles : angles equisitantly repartitioned along the arc

        """
        object_total_arc= self.func_arc(domain[1]) - self.func_arc(domain[0]) 
        angles=np.zeros([self.N])
        for i in range(0,self.N-1):
            local_func = lambda t : self.func_arc(t) - self.func_arc(angles[i]) - object_total_arc/self.N 
            res = scipy.optimize.root_scalar(local_func,bracket=[angles[i],2*np.pi],xtol=1e-12)
            angles[i+1] = res.root
        return angles
    
    def two_dimensional_equidistant_points(self,include_normals=False):
        """
        Computes the boundary points of the domain enclosed by the arc p.
        This function assumes that the domain is centered around the origin.

        Parameters
        ----------
        include_normals : boolean, optional
            Computes and returns the exterior unit normals at the points X if set to True.
            The default is False.

        Returns
        -------
        X : numpy array
        returns a 2D numpy array containing points on the boundary.
        normals : numpy array
        returns a 2D numpy array containing the normals vectors.
        """
        try:
            angles=self.angles_gen_arc(self.fourier_coeff,self.domain)
        except:
            raise RuntimeError("No function associated with the arc. Call the method fit_fourier before calling two_dimensionam_equidistant_points.")
        X=np.zeros([self.N,2])
        X[:,0]= np.cos(angles)*self.func_p(angles) -np.sin(angles)*self.func_p_prime(angles)
        X[:,1]= np.sin(angles)*self.func_p(angles)+np.cos(angles)*self.func_p_prime(angles)
        self.X=X
        if include_normals:
            X1= np.roll(X,-1,axis=0)
            X2=np.roll(X,1,axis=0)
            first= X1-X
            second=X-X2
            f_perp = np.array([-first[:, 1], first[:, 0]]).T
            s_perp = np.array([-second[:, 1], second[:, 0]]).T
            normals= (f_perp+s_perp)/2
            norms= np.linalg.norm(normals,axis=1)
            normals= normals/ np.reshape(norms,(-1,1))
            return X, normals
        else:
            return X
        
    def two_dimensional_plotter(self, plot_type="scatter"):
        """
        Create a 2D plot based on the specified plot type.
    
        Parameters:
            plot_type (str, optional): Type of plot to create ("scatter" or "line"). Default is "scatter".
        """
        Z = np.zeros([self.N + 1, 2])
        Z[:-1, :] = self.X
        Z[-1, :] = self.X[0, :]
    
        if plot_type == "scatter":
            plt.scatter(self.X[:, 0], self.X[:, 1])
        elif plot_type == "line":
            plt.plot(Z[:, 0], Z[:, 1])
        else:
            raise ValueError("Invalid plot_type. Choose either 'scatter' or 'line'.")
    
        plt.title(f"{plot_type.capitalize()} Plot")
        plt.xlabel("X-axis")
        plt.ylabel("Y-axis")
        plt.show()
    
class point_selector_dimension_three():
    
    def __init__(self):
        pass 
    
    def fit_spherical_funcs(self,r,r_theta,r_phi,N,inter_t=(0,2*np.pi),inter_s=(0,np.pi)):
        """
        Fit spherical functions for a 3D surface.
      
        Parameters:
            r (callable): The spherical function defining the radial distance.
            r_theta (callable): The derivative of the spherical function with respect to theta.
            r_phi (callable): The derivative of the spherical function with respect to phi.
            N (int): The maximal number of points on the largest band of the domain.
            inter_t (tuple): Tuple representing the interval for parameter theta (start_t, end_t).
            inter_s (tuple): Tuple representing the interval for parameter phi (start_s, end_s).

      
        Returns:
            None
      
        Description:
            This method fits spherical functions for a 3D surface based on the provided spherical coordinates.
            The parameters 'r', 'r_theta', and 'r_phi' are callables representing the spherical function,
            its derivative with respect to theta, and its derivative with respect to phi, respectively.
            'N' signifies the maximal number of points on each band of the domain.
      
            The fitted functions are stored as attributes: 'self.r', 'self.r_theta', 'self.r_phi'.
        """
        self.r = r
        self.r_theta = r_theta
        self.r_phi= r_phi
        self.N = N
        self.inter_t=inter_t
        self.inter_s = inter_s
        self.tracker="spherical"
        if N <=0:
            raise ValueError("The variable N needs to be strictly positive.")
        elif inter_t[1] <= inter_t[0]:
            raise ValueError("Invalid interval for the angle t.")
        elif inter_s[1] <= inter_s[0]:
            raise ValueError("Invalid interval for the angle s.")
        elif inter_t[1] <= inter_t[0]:
            raise ValueError("Invalid interval for the first angle.")
        
    def fit_parametric_funcs(self,f,f_theta,f_phi,N,inter_t=(0,2*np.pi),inter_s=(0,np.pi)):
        """
        Fit parametric surface functions for a 3D surface.
      
        Parameters:
            f (callable): Parametric function representing the surface.
            f_theta (callable): Derivative of f with respect to theta.
            f_phi (callable): Derivative of f with respect to phi.
            N (int):The maximal number of points on the largest band of the domain.
            inter_t (tuple): Tuple representing the interval for parameter theta (start_t, end_t).
            inter_s (tuple): Tuple representing the interval for parameter phi (start_s, end_s).
      
        Returns:
            None
      
        Description:
            This method fits parametric surface functions for a 3D surface.
            The parameters 'f', 'f_theta', and 'f_phi' are callables representing the parametric surface function,
            its derivative with respect to theta, and its derivative with respect to phi, respectively.
            The maximal number of points on the largest band of the domain.
      
            The fitted functions are stored as attributes: 'self.f', 'self.f_theta', 'self.f_phi'.
        """
        self.f = f
        self.f_theta = f_theta
        self.f_phi= f_phi
        self.N = N
        self.inter_t=inter_t
        self.inter_s = inter_s
        self.tracker="parametric"
        if N <=0:
            raise ValueError("The variable N needs to be strictly positive.")
        elif inter_t[1] <= inter_t[0]:
            raise ValueError("Invalid interval for the angle t.")
        elif inter_s[1] <= inter_s[0]:
            raise ValueError("Invalid interval for the angle s.")
        elif inter_t[1] <= inter_t[0]:
            raise ValueError("Invalid interval for the first angle.")
        
    def fit_harmonics(self,params_lm,L,N):
        """
        Fit spherical harmonics to the given parameters.
    
        Parameters:
            params_lm (numpy.ndarray): Array containing the spherical harmonic coefficients.
            L (int): The maximum degree of the spherical harmonics.
    
        Returns:
            None
    
        Description:
            This method fits spherical harmonics to the given parameters. The spherical harmonic
            coefficients are provided in the 'params_lm' array, and 'L' represents the maximum
            degree of the spherical harmonics.
    
            The fitted harmonics are stored as attributes: 'self.params_lm' and 'self.L'.
        """
        self.params_lm = params_lm
        self.L =L
        self.N = N
        self.inter_t=(0,2*np.pi)
        self.inter_s=(0,np.pi)
        self.tracker="harmonical"
        if N <=0:
            raise ValueError("The variable N needs to be strictly positive.")
        elif len(params_lm) != (L+1)**2:
            raise ValueError(f"The length of the array of spherical coefficients (here {len(params_lm)}) needs to be equal to (L+1)**2 (here {(L+1)**2}).")
    
    
    def main_point_generator(self):
        """
        Main method for generating points and normals based on the specified tracker type.
    
        This method checks the 'tracker' attribute and calls the appropriate point_generator function based on the tracker type:
        - If tracker is "parametric," it uses the provided parametric functions and interpolations.
        - If tracker is "spherical," it uses the spherical functions and interpolations.
        - If tracker is "harmonical," it uses the harmonical functions and interpolations.
    
        Returns:
        - X: Generated points.
        - normals: Generated normals.
    
        Note: Ensure that the required functions (point_generator, f, f_theta, f_phi, r, r_theta, r_phi, creator_r, creator_r_theta, creator_r_phi) and interpolations (inter_t, inter_s) are defined in the class.
    
        Raises:
        - RuntimeError: Raised if no function associated with the surface. Call any of the fit_polar_funcs, fit_spherical_funcs, or fit_harmonics methods before main_initialization.
       """
        try:
            self.tracker
        except:
            raise RuntimeError("No function associated with the surface. Call any of the fit_polar_funcs, fit_spherical_funcs, or fit_harmonics methods before main_initialization.")
        if self.tracker == "parametric":
            X, normals = self.point_generator(self.f,self.f_theta,self.f_phi,self.inter_t,self.inter_s)
        elif self.tracker == "spherical":
            f = lambda t, s : self.r(t,s) * np.array([ np.cos(t)*np.sin(s) , np.sin(t)* np.sin(s),np.cos(s) ])
            f_t = lambda t, s : self.r(t,s)* np.array([-np.sin(t)*np.sin(s) , np.cos(t)* np.sin(s),0 ]) + self.r_theta(t,s) *np.array([ np.cos(t)*np.sin(s) , np.sin(t)* np.sin(s),np.cos(s) ])
            f_s = lambda t,s : self.r(t,s) *    np.array([np.cos(t)*np.cos(s) , np.sin(t)* np.cos(s),-np.sin(s)]) + self.r_phi(t,s) * self.r(t,s) * np.array([ np.cos(t)*np.sin(s) , np.sin(t)* np.sin(s),np.cos(s) ])
            X, normals = self.point_generator(f,f_t,f_s,self.inter_t,self.inter_s)
        elif self.tracker == "harmonical":
            f = lambda t, s : self.creator_r(t,s) * np.array([ np.cos(t)*np.sin(s) , np.sin(t)* np.sin(s),np.cos(s) ])
            f_t = lambda t, s : self.creator_r(t,s)* np.array([-np.sin(t)*np.sin(s) , np.cos(t)* np.sin(s),0 ]) + self.creator_r_theta(t,s) *np.array([ np.cos(t)*np.sin(s) , np.sin(t)* np.sin(s),np.cos(s) ])
            f_s = lambda t,s : self.creator_r(t,s) *    np.array([np.cos(t)*np.cos(s) , np.sin(t)* np.cos(s),-np.sin(s)]) + self.creator_r_phi(t,s) * self.creator_r(t,s) * np.array([ np.cos(t)*np.sin(s) , np.sin(t)* np.sin(s),np.cos(s) ])
            X, normals = self.point_generator(f,f_t,f_s,self.inter_t,self.inter_s)
        
        return X, normals


    def point_generator(self,f,f_theta,f_phi,inter_t,inter_s):
        """
        Generate 3D points on a parametric surface defined by functions f, f_theta, and f_phi.
    
        Parameters:
            f (callable): Parametric function representing the surface.
            f_theta (callable): Derivative of f with respect to theta.
            f_phi (callable): Derivative of f with respect to phi.
            inter_t (tuple): Tuple representing the interval for parameter theta (start_t, end_t).
            inter_s (tuple): Tuple representing the interval for parameter phi (start_s, end_s).
    
        Returns:
            tuple: A tuple containing two numpy arrays:
                - First array: 3D points representing the generated surface.
                - Second array: Normalized normals at each generated point.
    
        Description:
            This method generates equidistant points on a 3D surface defined by parametric functions.
            The parameters 'f', 'f_theta', and 'f_phi' are callables representing the surface function,
            its derivative with respect to theta, and its derivative with respect to phi, respectively.
            'inter_t' and 'inter_s' are tuples specifying the intervals for theta and phi, respectively.
    
            The method returns a tuple containing two arrays:
            - The first array represents the 3D points on the surface.
            - The second array contains the normalized normals corresponding to each generated point.
         """
        start_t= inter_t[0]
        end_t= inter_t[1]
        start_s= inter_s[0]
        end_s= inter_s[1]
        curv_0 = lambda t :np.linalg.norm( f_phi(start_t,t))
        arc_length = scipy.integrate.quad(curv_0,start_s,end_s)[0]/self.N
        points_s = np.zeros([self.N,3])
        list_s = np.zeros([self.N])
        points_s[0,]=f(start_t,start_s)
        list_s[0]=start_s
        spread_s= end_s - start_s
        spread_t = end_t - start_t
        for i in range(0,self.N-1):
            local_func = lambda s : np.linalg.norm( f(start_t,s) - points_s[i,:]) - arc_length
            res = scipy.optimize.root_scalar(local_func,x0= list_s[i]+ spread_s/self.N-0.01, x1=spread_s*(i+1)/self.N)
            points_s[i+1,:] = f(start_t,res.root)
            list_s[i+1]=res.root
        complete_pts=np.zeros([0,3])
        complete_normals_t = np.zeros([0,3])
        complete_normals_s = np.zeros([0,3])
        for s in list_s:
            loc_curv = lambda t :np.linalg.norm(f_theta(t,s))
            arc_length_loc = scipy.integrate.quad(loc_curv,start_t,end_t)[0]
            if arc_length_loc !=0:
                number_pts= int(arc_length_loc/arc_length)+1
                local_pts=np.zeros([number_pts,3])
                local_pts[0,:]= f(start_t,s)
                local_norm_t=np.zeros([number_pts,3])
                local_norm_s=np.zeros([number_pts,3])
                local_list=np.zeros([number_pts])
                local_list[0]=0
                local_norm_t[0]= f_theta(start_t,s)
                local_norm_s[0]=  f_phi(start_t,s)
                for j in range(0,number_pts-1):
                    local_func = lambda t : np.linalg.norm( f(t,s) - local_pts[j,:]) - arc_length_loc/number_pts
                    res = scipy.optimize.root_scalar(local_func,x0= local_list[j] +spread_t /number_pts-0.01,x1= local_list[j] +spread_t /number_pts+0.01)
                    local_pts[j+1,:]=f(res.root,s)
                    local_norm_t[j+1,:]= f_theta(res.root,s)
                    local_norm_s[j+1,:]= f_phi(res.root,s)
                    local_list[j+1]=res.root
                complete_pts = np.concatenate((complete_pts,local_pts))
                complete_normals_t = np.concatenate((complete_normals_t,local_norm_t))
                complete_normals_s = np.concatenate((complete_normals_s,local_norm_s))
        normals= np.cross(complete_normals_s, complete_normals_t)
        norms= np.linalg.norm(normals,axis=1)
        normals = normals/np.reshape(norms,(-1,1))
        self.X = complete_pts   
        self.normalized_normals = normals
        return complete_pts,normals
    def shape_plotter(self,ann=False):
        
        ax = plt.axes(projection ='3d')
        '''
        Z[:self.N,:]=self.X
        Z[self.N,0]=self.X[0,0]
        Z[self.N,1]=self.X[0,1]
        Z[self.N,2]=self.X[0,2]
        '''
        ax.scatter3D(self.X[:,0], self.X[:,1], self.X[:,2],color= 'green')
        if ann==True:
            
            for i in range(self.N):
                ax.text(self.X[i,0],self.X[i,1],self.X[i,2],i)
    
    
    
    