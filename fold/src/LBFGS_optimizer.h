/*
===============================================================================
   Implementation of the LBFGS optimization algorithm in C/C++   
   
   Some code is referenced from Rosetta (https://www.rosettacommons.org/software)
           
   Please report bugs and questions to robpearc@umich.edu
===============================================================================
*/

#ifndef LBFGS_OPTIMIZER
#define LBFGS_OPTIMIZER

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <memory>

#include "Operations.h"
#include "CommonParameters.h"
#include "Energy.h"

extern Energy energyFunction;

using namespace std;

bool fractional_converge_test( double Fnew, double Fold )
{	
    return ( 2.0f * std::abs( Fnew - Fold ) <=
             1.0e-4 * ( std::abs( Fnew ) + std::abs( Fold ) + 1.0e-10 ));	
}

bool absolute_converge_test( double Fnew, double Fold )
{
    return ( std::abs( Fnew - Fold ) <= 1.0e-2 );
}

///////////////////////////////////////////////////////////
// Class for data stored at each LBFGS iteration
///////////////////////////////////////////////////////////
class lbfgs_iteration_data {
public:
    double alpha;
    double *s;
    double *y;
    double ys; // = dot(y,s)

    ~lbfgs_iteration_data(){
        delete [] s;
        delete [] y;
    }
};

///////////////////////////////////////////////////////////
// Base class for reference-counted polymorphic classes
///////////////////////////////////////////////////////////
class ReferenceCount1
{
public: // Creation

    /// @brief Default constructor
    inline
    ReferenceCount1()
    {}

    inline
    virtual
    ~ReferenceCount1()
    {}
}; // ReferenceCount1

///////////////////////////////////////////////////////////
// Functor that stores current position and
// search direction to provide a double->double mapping
// for univariate minimization
///////////////////////////////////////////////////////////
class func_1d {
private:
    int _size;
    double *_starting_point;
    double *_search_direction;
    double *_eval_point;
public:
    double *_dE_dvars;
private:
    int _eval_count;
    int _deriv_count;
public:
    func_1d( int size, double *start, double *dir  ) :
        _size(size), 
        _starting_point( start ),
        _search_direction( dir ),
        _eval_point ( new1DArr( _size ) ),
        _dE_dvars ( new1DArr( _size ) ),
        _eval_count( 0 ),
        _deriv_count( 0 )
    {
    };
	
    virtual ~func_1d();

    double operator()( double displacement, point3f *decstr )
    {
        _eval_count++;
        for ( uint i = 0 ; i <_size ; ++i ) 
        {
            _eval_point[i] = _starting_point[i] + ( displacement * _search_direction[i] );
        }
        return energyFunction.calcrmsdenergy( decstr , _eval_point );
    };

    void reset_eval_count() { _eval_count = 0; };
    int get_eval_count() { return _eval_count; };
    int get_deriv_count() { return _deriv_count; };
    /// @brief Error condition wherein the computed gradient does not match the actual gradient;
    /// invokes the Multifunc::dump( vars, vars2 ) method.
    void dump( double displacement );
};

func_1d::~func_1d() = default;

void func_1d::dump( double displacement ) {
	for ( int i = 0 ; i < _size ; i++ ) {
		_eval_point[i] = _starting_point[i] + ( displacement * _search_direction[i] );
	}
	//return _func.dump( _starting_point, _eval_point );
}

///////////////////////////////////////////////////////////
// base class / interface for line minimizers
///////////////////////////////////////////////////////////
class LineMinimizationAlgorithm : public ReferenceCount1
{
public:
    ~LineMinimizationAlgorithm() override;
    LineMinimizationAlgorithm( int dimension ) :
        _last_accepted_step( 1.0 ), _func_to_beat( 0.0 ),
        _deriv_sum( 0.0 ), _num_linemin_calls(0), _tolerance( 0.1 ),
        _nonmonotone( false ), _silent( true ) {};
    virtual double operator()( point3f *, int, double *, double * ){ return 0.0; };
    virtual bool provide_stored_derivatives(){ return false; };
    bool nonmonotone() { return _nonmonotone; };

    bool silent() { return _silent; };
    void silent(bool s_in) { _silent=s_in; };

    double _last_accepted_step;
    double _func_to_beat;
    double _deriv_sum;
    int _num_linemin_calls;

protected:
    double const _tolerance;
    bool _nonmonotone;
    bool _silent;
};

/// @details Auto-generated virtual destructor
LineMinimizationAlgorithm::~LineMinimizationAlgorithm() = default;

///////////////////////////////////////////////////////////
// concrete line minimizer - Armijo's method
///////////////////////////////////////////////////////////
class ArmijoLineMinimization : public LineMinimizationAlgorithm
{
public:
    ArmijoLineMinimization( bool nonmonotone, int dim, double max_step_limit ) :
        LineMinimizationAlgorithm( dim ),
        _num_calls( 0 ),
        max_step_limit_( max_step_limit ) { _nonmonotone = nonmonotone; };

    bool provide_stored_derivatives() override{ return false; };
    double operator()( point3f *decstr, int dim, double * curr_pos, double * curr_dir ) override;
    double Armijo( double init_step, point3f *decstr, func_1d & func_eval );

    int _num_calls;
    double max_step_limit_;
};


double ArmijoLineMinimization::operator()( point3f *decstr, int dim, double *current_position, double *search_direction )
{
    double const FACTOR( 0.5 );
    int const problem_size( dim );

    _num_linemin_calls++;

    //cout<<"ArmijoLineMinimization::operator()"<<endl;
    // Construct the one-dimensional projection of the function
    func_1d this_line_func( problem_size, current_position, search_direction );

    // Early termination for derivatives (search_direction magnitudes) near zero
    // Please note that the search_direction vector is not normalized
    double derivmax = 0.0;
    for ( int i = 0 ; i < problem_size; i++ ) 
    {
        if ( std::abs( search_direction[ i ] ) >
             std::abs( derivmax ) ) derivmax = search_direction[ i ];
    }

    // if ( runlevel >= gush) std::cout << "derivmax," << SS( derivmax ) << std::endl;
    if ( std::abs(derivmax) < .0001 )
    {
        double final_value = this_line_func( 0.0, decstr ); // deriv = 0, return value
        return final_value;
    }

    //initial trial stepsize
    double init_step( _last_accepted_step / FACTOR );
    if ( init_step > max_step_limit_ ) init_step = max_step_limit_;

    double final_value = Armijo( init_step, decstr, this_line_func );

    for ( int j = 0; j < problem_size; j++ )
    {
        search_direction[j] *= _last_accepted_step;
        current_position[j] += search_direction[j];
    }

    // std::cout << "Linemin used " << this_line_func.get_eval_count() <<
    //  " function calls and returns " << final_value << " on step of " << _last_accepted_step << std::endl;

    return final_value;
}

double ArmijoLineMinimization::Armijo( double init_step, point3f *decstr, func_1d & func_eval ) 
{
    double const FACTOR( 0.5 );
    double const SIGMA( 0.1 );
    double const SIGMA2( 0.8 );
    static double const MINSTEP = 1e-9;

    //std::cout << "func_to_beat is " << _func_to_beat << std::endl;

    double func_value = func_eval( init_step, decstr );
    _num_calls++;
    _last_accepted_step = init_step;

    if ( func_value < _func_to_beat + init_step * SIGMA2 * _deriv_sum )
    {
        double test_step = init_step/FACTOR;
        double test_func_value = func_eval( test_step, decstr );
        _num_calls++;
        if ( test_func_value < func_value )
	{
            _last_accepted_step = test_step;
            return test_func_value;
        }
        return func_value;
    }

    //cout << "init_step " << init_step << endl;
    double far_step = init_step;
    while ( func_value > _func_to_beat + init_step*SIGMA*_deriv_sum ) 
    {
        // Abort if function value is unlikely to improve.
        if ( ( init_step <= 1e-5 * far_step ) ||
             (init_step < MINSTEP && func_value >= _func_to_beat) ) 
        {
             
            double test_step = ( func_value - _func_to_beat ) / init_step;
            if ( !_silent ) 
            {
                cout << "Inaccurate G! step= " << ( init_step ) << " Deriv= " <<
                        ( _deriv_sum ) << " Finite Diff= " << ( test_step ) << endl;
            }
            func_eval.dump( init_step );
            _last_accepted_step = 0.0;
            return _func_to_beat;
        }

        init_step *= FACTOR*FACTOR;  // faster step decrease
        func_value = func_eval( init_step, decstr );
        _num_calls++;
    }

    _last_accepted_step = init_step;

    if ( init_step < 0.0 ) 
    {
        cout << "Forced to do parabolic fit!" << std::endl;
        // Parabola interpolate between 0 and init_step for refinement
        double test_step = -_deriv_sum*init_step*init_step/
                           (2*(func_value - _func_to_beat - init_step * _deriv_sum));
        if ( test_step > 1e-3*far_step && test_step < far_step ) 
        {
            double test_func_value = func_eval( test_step, decstr );
            _num_calls++;
            if ( test_func_value < func_value ) 
            {
                _last_accepted_step = test_step;
                func_value = test_func_value;
            }
        }
    }

    return func_value;
}


///////////////////////////////////////////////////////////
// Forward
///////////////////////////////////////////////////////////
class LineMinimizationAlgorithm;
typedef std::shared_ptr< LineMinimizationAlgorithm > LineMinimizationAlgorithmOP;


///////////////////////////////////////////////////////////
///Simple low-level minimizer class
///////////////////////////////////////////////////////////
class Minimizer 
{
public:
    Minimizer();
    double run( double *vars_inout, point3f *decstr, int dim, int ITMAX, string ConvergeTest );

    ~Minimizer()
    {
        delete [] XP; XP = NULL;
        delete [] Xtemp; Xtemp = NULL;		
        delete [] G; G = NULL;
        delete [] GP; GP = NULL;
        delete [] Gtemp; Gtemp = NULL;
        delete [] D; D = NULL;
        delete [] W; W = NULL;
        delete [] pf; pf = NULL;
        delete [] lm;
    };


private:
	void lbfgs(
            double *X,
        point3f *decstr,
        int dim,
        double &FRET,
        LineMinimizationAlgorithmOP line_min,
        int const ITMAX,
        string ConvergeTest,
        bool w_rescore = false /*= false*/);

        
private:
	double *XP;					//current position
	double *Xtemp;		
	double *G;					//gradient vectors Gk
	double *GP;					//current gradient vectors Gk+1
	double *Gtemp;
	double *D;					//direction
	double *W;					//?
	double *pf;
	lbfgs_iteration_data *lm;	
}; 


// set the function and the options
Minimizer::Minimizer()
{
}

double Minimizer::run( double *vars_inout, point3f *decstr, int dim, int ITMAX, string ConvergeTest )
{
    double * vars( vars_inout );
    double end_func;
    double const ARMMAXSTEP( 20.0 );
    //int const ITMAX( 200 );
    
    LineMinimizationAlgorithmOP armijo_line_search( new ArmijoLineMinimization( true, dim, ARMMAXSTEP ) );
    armijo_line_search->silent( true );
    lbfgs( vars, decstr, dim, end_func, armijo_line_search, ITMAX, ConvergeTest, false );
    vars_inout = vars;

    return energyFunction.calcrmsdenergy( decstr , vars );
}


////////////////////////////////////////////////////////////////////////
// *      Limited memory BFGS (L-BFGS).
// *
// * Copyright (c) 1990, Jorge Nocedal
// * Copyright (c) 2007-2010 Naoaki Okazaki
// * All rights reserved.
// *
// * Permission is hereby granted, free of charge, to any person obtaining a copy
// * of this software and associated documentation files (the "Software"), to deal
// * in the Software without restriction, including without limitation the rights
// * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// * copies of the Software, and to permit persons to whom the Software is
// * furnished to do so, subject to the following conditions:
// *
// * The above copyright notice and this permission notice shall be included in
// * all copies or substantial portions of the Software.
// *
// * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// * THE SOFTWARE.
////////////////////////////////////////////////////////////////////////

// This library is a C port of the FORTRAN implementation of Limited-memory
// Broyden-Fletcher-Goldfarb-Shanno (L-BFGS) method written by Jorge Nocedal.
// The original FORTRAN source code is available at:
// http://www.ece.northwestern.edu/~nocedal/lbfgs.html
//
// The L-BFGS algorithm is described in:
//     - Jorge Nocedal.
//       Updating Quasi-Newton Matrices with Limited Storage.
//       <i>Mathematics of Computation</i>, Vol. 35, No. 151, pp. 773--782, 1980.
//     - Dong C. Liu and Jorge Nocedal.
//       On the limited memory BFGS method for large scale optimization.
//       <i>Mathematical Programming</i> B, Vol. 45, No. 3, pp. 503-528, 1989.
void Minimizer::lbfgs( double *X, point3f *decstr, int dim, double &FRET,
                       LineMinimizationAlgorithmOP line_min, int const ITMAX,
                       string ConvergeTest, bool w_rescore /*= false*/) 
{	
    int const N( dim );
    static int M( 256 );
    int const PAST( line_min->nonmonotone() ? 3 : 1 );
    double const EPS( 1.E-5 );
    double gmax_cutoff_for_convergence = 1e-8;
    bool silent = true;

    int K = 1; // number of func evaluations
    XP = new1DArr( N );          //current position
    Xtemp = new1DArr( N );
    G = new1DArr( N );							//gradient vectors Gk
    GP = new1DArr( N );							//current gradient vectors Gk+1
    Gtemp = new1DArr( N );
    D = new1DArr( N );							//direction
    W = new1DArr( N );

    // Allocate & initialize limited memory storage
    int CURPOS = 1; 							// pointer to current location in lm
    lm = new lbfgs_iteration_data[M];
    for ( int i=0; i < M; i++ )
    {
        lm[i].alpha = 0;
        lm[i].ys = 0;
        lm[i].s = new1DArr(N);
        lm[i].y = new1DArr(N);
    }

    // Allocate space for storing previous values of the objective function
    pf = new1DArr( PAST );

    // Evaluate the function value and its gradient
    int func_memory_filled( 1 );
    double prior_func_value = energyFunction.calcrmsdenergy( decstr, X, G, true );//func_( X );
    pf[0] = FRET = prior_func_value;

    // Compute the direction
    // assume the initial hessian matrix H_0 as the identity matrix.
    double invdnorm = 0.0;
    for ( int i = 0; i < N; i++ ) {
        D[i] = -G[i];
    }

    if ( line_min->nonmonotone() ) line_min->_last_accepted_step = 0.005;
	
    bool last_step_good = true;
    for ( int ITER = 1; ITER <= ITMAX; ++ITER ) 
    {
        if(ITER%100==0) cout << "Step: " << ITER << " ";
        // Store the current position and gradient vectors
        if ( last_step_good ) {
            for(int i=0; i<N; i++){
                XP[i] = X[i];
                GP[i] = G[i];
            }
        }

        // line min
        line_min->_deriv_sum = 0.0;
        double Gmax = 0.0;
        for ( int i = 0; i < N; i++) {
            line_min->_deriv_sum += D[i] * G[i];
            if ( std::abs( G[i] ) > Gmax ) {
                Gmax = std::abs( G[i] );
            }
        }

        line_min->_func_to_beat = pf[ 0 ];
        for ( int i = 1 ; i < func_memory_filled ; i++ ) {
            if ( line_min->_func_to_beat < pf[ i ] ) {
                line_min->_func_to_beat = pf[ i ];
            }
        }
        if(ITER%100==0) cout << "Energy: "<<line_min->_func_to_beat << endl;
        // check 1: if derivative is positive, flip signs of positive components
        if ( line_min->_deriv_sum > -EPS ) {
            line_min->_deriv_sum = 0.0;
            for ( int i = 0; i < N; i++ ) {
                if ( D[i]*G[i] >= 0 ) D[i] = -D[i];

                line_min->_deriv_sum += D[i]*G[i];
                if ( std::abs( G[i] ) > Gmax ) {
                    Gmax = std::abs( G[i] );
                }
            }
        }

        // check 2: if derivative still positive, reset Hessian
        if ( line_min->_deriv_sum > -EPS ) {
            line_min->_deriv_sum = 0.0;
            for ( int i = 0; i < N; i++ ) {
                D[i] = -G[i];
                line_min->_deriv_sum += D[i]*G[i];
            }

            if ( sqrt( -line_min->_deriv_sum ) > 1e-6 ) {
                func_memory_filled = 1;
                line_min->_func_to_beat = pf[0] = prior_func_value = FRET;
            }
        }

        // X is returned as new pt, and D is returned as the change
        FRET = (*line_min)( decstr, N, X, D );

        bool covergeTest = false;
        if ( ConvergeTest == "fractional"){
            covergeTest = fractional_converge_test( FRET, prior_func_value );
        }
        else{
            covergeTest = absolute_converge_test ( FRET, prior_func_value );
        }
        if ( covergeTest || line_min->_last_accepted_step == 0 ) {
            if ( Gmax <= gmax_cutoff_for_convergence ) {
                return;
            } 
	    else {
                if ( line_min->_last_accepted_step == 0 ) { // failed line search
                    // Reset Hessian
                    CURPOS = 1;
                    K = 1;

                    // reset line minimizer
                    line_min->_deriv_sum = 0.0;
                    for ( int i = 0; i < N; i++ ) {
                        D[i] = -G[i];
                        line_min->_deriv_sum += D[i]*G[i];
                    }
  
                    if ( line_min->_deriv_sum <= -EPS ) {
                        invdnorm = 1.0/sqrt( -line_min->_deriv_sum );

                        // delete prior function memory
                        line_min->_last_accepted_step = invdnorm; 
                        // reset initial step size (in case default of '1' is too large)
                        func_memory_filled = 1;
                        prior_func_value = FRET;
                        pf[0] = prior_func_value;

                        //line search in the direction of the gradient
                        FRET = (*line_min)( decstr, N, X, D );
                    } 
                    else {
                        return;
                    }

                    // if the line minimzer fails again, abort
                    if ( line_min->_last_accepted_step == 0 ) {
                        if ( !silent ) {
                            cout << "Line search failed even after resetting Hessian; aborting at iter#" << ITER << endl;
                        }
                        return;
                    }
                }
            } 
        }

        prior_func_value = FRET;

        // Update memory of function calls
        if ( func_memory_filled < PAST ) {
            func_memory_filled++;
        }
        else {
            for ( int i = 0 ; i < PAST-1 ; i++ ) {
                pf[ i ] = pf[ i + 1 ];
            }
        }
        pf[ func_memory_filled - 1] = prior_func_value;

        for ( int i = 0; i < N; i++ ) {
            GP[i] = G[i];
        }

        energyFunction.calcrmsdenergy( decstr, X, G, true );

        line_min->_deriv_sum = 0.0;
        double deriv_new = 0.0;
        for ( int i = 0; i < N; i++) {
            line_min->_deriv_sum += D[i]*GP[i];
            deriv_new += D[i]*G[i];
        }

        // LBFGS updates
        //
        // Update vectors s and y: 
        //   s_{k+1} = x_{k+1} - x_{k} = \step * d_{k}.
        //   y_{k+1} = g_{k+1} - g_{k}.
        // Compute scalars ys and yy:
        //   ys = y^t \cdot s = 1 / \rho.
        //   yy = y^t \cdot y.
        // Notice that yy is used for scaling the hessian matrix H_0 (Cholesky factor).
        double ys=0, yy=0;
        for ( int i = 0; i < N; i++ ) {
            Xtemp[i] = X[i] - XP[i];
            Gtemp[i] = G[i] - GP[i];
            ys += Gtemp[i]*Xtemp[i];
            yy += Gtemp[i]*Gtemp[i];
        }

        if ( std::fabs( ys ) < 1e-6 ) {
            last_step_good = false;
        } else {
            last_step_good = true;
            for ( int i = 0; i < N; ++i ) {
                lm[CURPOS-1].s[i] = Xtemp[i];
                lm[CURPOS-1].y[i] = Gtemp[i];
            }
            lm[CURPOS-1].ys = ys;

            // increment
            K++;
            CURPOS++;
            if ( CURPOS > M ) CURPOS = 1;
        }

        // Recursive formula to compute dir = -(H \cdot g).
        //   This is described in page 779 of:
        //   Jorge Nocedal.
        //   Updating Quasi-Newton Matrices with Limited Storage.
        //   Mathematics of Computation, Vol. 35, No. 151,
        //   pp. 773--782, 1980.
        int bound = std::min( M, K-1 );

        // Compute the negative of gradients
        for ( int i = 0; i < N; i++ ) {
            D[i] = -G[i];
        }

        //backward
        int j = CURPOS-1;
        for ( int pts=0; pts<bound; pts++ ) {
            j--;
            if ( j<0 ) j=M-1; // wrap around

            // \alpha_{j} = \rho_{j} s^{t}_{j} \cdot q_{k+1}
            lm[j].alpha = 0;
            for ( int i = 0; i < N; i++ ) {
                lm[j].alpha += lm[j].s[i] * D[i];
            }
            lm[j].alpha /= lm[j].ys;

            // q_{i} = q_{i+1} - \alpha_{i} y_{i}
            for ( int i = 0; i < N; i++ ) {
                D[i] += -lm[j].alpha * lm[j].y[i];
            }
        }
 
        //forward
        for ( int pts=0; pts<bound; ++pts ) {
            //if (std::fabs(lm[j].ys) < 1e-6) continue;

            // \beta_{j} = \rho_{j} y^t_{j} \cdot \gamma_{i}
            double beta=0.0;
            for ( int i = 0; i < N; i++ ) {
                beta += lm[j].y[i] * D[i];
            }
            beta /= lm[j].ys;

            // \gamma_{i+1} = \gamma_{i} + (\alpha_{j} - \beta_{j}) s_{j}
            for ( int i = 0; i < N; i++ ) {
                D[i] += (lm[j].alpha - beta) * lm[j].s[i];
            }

            j++;
            if ( j>M-1 ) j=0; // wrap around
        }
    }

    if ( !silent ) {
        cout << "Warning! LBFGS MAX CYCLES " << ITMAX << " EXCEEDED, BUT FUNC NOT CONVERGED!" << endl;
    }
	
    return;
}

#endif
