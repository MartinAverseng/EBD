classdef H0Kernel < Kernel
    properties
        k; 
        C;
    end
    methods
        function[this] = H0Kernel(kk,CC)
            if nargin == 0
                kk = 1;
            end
            if nargin <= 1
                CC = 1;                                
            end
            this.k = kk;
            this.C = CC;
            modelKern = this.C*(J0Kernel(kk) + 1i*Y0Kernel(kk));
            this.func = modelKern.func;
            this.der = modelKern.der;
            this.scalFunc = modelKern.scalFunc;
        end
        function[out] = dilatation(this,lambda)
            out = H0Kernel(this.k*lambda,this.C);
        end
        function[out] = mtimes(this,mu)
            if isa(this,'Kernel')
                assert(and(isa(mu,'double'),isscalar(mu)));
                out = H0Kernel(this.k,this.C*mu);
            else
                out = mtimes(mu,this);
            end
        end
        function[rq] = radialQuadKernel(this,a,tol,varargin)
            modelKern = J0Kernel(this.k) + 1i*Y0Kernel(this.k);
            rqTemp = modelKern.radialQuadKernel(a,tol/abs(this.C),varargin{:});
            rq = this.C*rqTemp;
        end
    end
end
