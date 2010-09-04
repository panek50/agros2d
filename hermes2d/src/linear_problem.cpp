// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "common.h"
#include "limit_order.h"
#include "discrete_problem.h"
#include "linear_problem.h"
#include "weakform.h"
#include "solver.h"
#include "space.h"
#include "precalc.h"
#include "refmap.h"
#include "solution.h"
#include "integrals_h1.h"
#include "views/view.h"
#include "views/vector_view.h"
#include "tuple.h"
#include "norm.h"


LinearProblem::LinearProblem() : DiscreteProblem() {};
LinearProblem::LinearProblem(WeakForm* wf_) : DiscreteProblem(wf_) {};
LinearProblem::LinearProblem(WeakForm* wf_, Space* s_) : DiscreteProblem(wf_, s_) {};
LinearProblem::LinearProblem(WeakForm* wf_, Tuple<Space*> spaces_) : DiscreteProblem(wf_, spaces_) {};
LinearProblem::~LinearProblem() {};

void LinearProblem::assemble(Matrix* mat_ext, Vector* rhs_ext, bool rhsonly, bool is_complex)
{
  int ndof = this->get_num_dofs();
  if (ndof == 0) error("ndof == 0 in LinearProblem::assemble().");
  Vector* dir_ext = new AVector(ndof, is_complex);
  // The vector dir represents the contribution of the Dirichlet lift, 
  // and for linear problems it has to be subtracted from the right hand side.
  // The NULL stands for the initial coefficient vector that is not used.
  DiscreteProblem::assemble(NULL, mat_ext, dir_ext, rhs_ext, rhsonly);
  // FIXME: Do we really need to handle the real and complex cases separately?
  if (is_complex) for (int i=0; i < ndof; i++) rhs_ext->add(i, -dir_ext->get_cplx(i));
  else for (int i=0; i < ndof; i++) rhs_ext->add(i, -dir_ext->get(i));
  delete dir_ext;
}


