#include <sundials/sundials_linearsolver.h>
#include "TChem.hpp"
#include "TChem_AerosolChemistry_CVODE_K_BatchedGMRES.hpp"

#define BGMR_CONTENT(S) (static_cast<TChem::SUNLinearSolverContent_BatchedGMRES*>((S)->content))

namespace TChem{

SUNLinearSolver SUNLinSol_BatchedGMRES(N_Vector y, 
                                       Impl::AerosolChemistryRHS<SUNLinearSolverContent_BatchedGMRES::problem_type> rhs_object,
                                       ordinal_type system_size, 
                                       ordinal_type n_systems,
                                       int pretype,
                                       ordinal_type gmres_max_iter,
                                       SUNLinearSolverContent_BatchedGMRES::policy_type policy,
                                       SUNContext sunctx)
{
  /* check that the supplied N_Vector supports all requisite operations */
  if (y == nullptr || y->ops->nvgetdevicearraypointer == nullptr || y->ops->nvgetlength == nullptr) {
    return nullptr;
  }
  // ensure the problem has the right stucture
  if (N_VGetLength(y) != (sunindextype)n_systems * system_size) { 
    return nullptr; 
  }

  /* check for legal pretype; if illegal use defaults */
  if ((pretype != SUN_PREC_NONE) && (pretype != SUN_PREC_RIGHT))
  {
    pretype = SUN_PREC_NONE;
  }

  if (gmres_max_iter <= 0) { 
    gmres_max_iter = BGMR_MAXITER_DEFAULT; 
  }

  // Create content
  SUNLinearSolverContent_BatchedGMRES* content;
  content = new SUNLinearSolverContent_BatchedGMRES(rhs_object, 
    system_size, n_systems, gmres_max_iter, pretype, policy);


  /* Create linear solver */
  SUNLinearSolver S = SUNLinSolNewEmpty(sunctx);
  if (S == nullptr) {
    delete content;
    return nullptr;
  }

  /* Attach operations */
  S->ops->gettype           = SUNLinSolGetType_BGMR;
  S->ops->getid             = SUNLinSolGetID_BGMR;
  S->ops->setatimes         = SUNLinSolSetATimes_BGMR;              // using to set cvode_mem off AData
  S->ops->setpreconditioner = SUNLinSolSetPreconditioner_BGMR;      // set pointers for Psetup, PSolve, Pdata
  S->ops->setscalingvectors = SUNLinSolSetScalingVectors_BGMR;      // just uses s1
  S->ops->setzeroguess      = SUNLinSolSetZeroGuess_BGMR;           // no-op shim, return success
  S->ops->initialize        = SUNLinSolInitialize_BGMR;             // calls SUNLinearSolverContent_BatchedGMRES::init()
  S->ops->setup             = SUNLinSolSetup_BGMR;                  // setup the preconditioner if user psetup routine passed to CVODE
  S->ops->solve             = SUNLinSolSolve_BGMR;                  // calls SUNLinearSolverContent_BatchedGMRES::solve()
  S->ops->numiters          = SUNLinSolNumIters_BGMR;               // return SUNLinearSolverContent_BatchedGMRES::numiters
  S->ops->resnorm           = SUNLinSolResNorm_BGMR;                // return SUNLinearSolverContent_BatchedGMRES::resnorm
  S->ops->resid             = NULL;                                 // unused
  S->ops->lastflag          = SUNLinSolLastFlag_BGMR;               // return SUNLinearSolverContent_BatchedGMRES::last_flag
  S->ops->space             = NULL;                                 // unused
  S->ops->free              = SUNLinSolFree_BGMR;                   // wrapper for calling .free()

  /* Attach content */
  S->content = content;

  return (S);
}

SUNErrCode SUNLinSol_BatchedGMRESSetVerbose(SUNLinearSolver S, bool onoff)
{
  if (S == nullptr || S->content == nullptr) { return SUN_ERR_ARG_CORRUPT; }
  BGMR_CONTENT(S)->verbose = onoff;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSol_BatchedGMRESSetStateDump(SUNLinearSolver S, const char* filename,
                                              long solve_index, sunrealtype rtol, sunrealtype atol)
{
  if (S == nullptr || S->content == nullptr) { return SUN_ERR_ARG_CORRUPT; }
  auto* c = BGMR_CONTENT(S);
  c->dump_file  = (filename == nullptr) ? "" : filename;
  c->dump_solve = solve_index;
  c->dump_rtol  = rtol;
  c->dump_atol  = atol;
  c->dumped     = false;
  return SUN_SUCCESS;
}

// Operations 
// TODO: Add exception catching inside each operator that isnt a no-op
SUNLinearSolver_Type SUNLinSolGetType_BGMR([[maybe_unused]] SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_ITERATIVE);
}

SUNLinearSolver_ID SUNLinSolGetID_BGMR([[maybe_unused]] SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_CUSTOM);
}

SUNErrCode SUNLinSolSetATimes_BGMR(SUNLinearSolver S, void* ATData,
                                    SUNATimesFn ATimes)
{
  /* set function pointers to integrator-supplied ATimes routine
     and data, and return with success */
  (void)ATimes; // unused
  BGMR_CONTENT(S)->cvode_mem = ATData;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetScalingVectors_BGMR(SUNLinearSolver S, N_Vector s1,
                                            N_Vector s2)
{
  /* set N_Vector pointers to integrator-supplied scaling vectors,
     and return with success */
  BGMR_CONTENT(S)->s1 = s1;
  //BGMR_CONTENT(S)->s2 = s2;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetZeroGuess_BGMR([[maybe_unused]] SUNLinearSolver S, [[maybe_unused]] sunbooleantype onff)
{
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolInitialize_BGMR(SUNLinearSolver S)
{
  auto* content = BGMR_CONTENT(S);
  return content->init();
}

SUNErrCode SUNLinSolSetPreconditioner_BGMR(SUNLinearSolver S, void* PData,
                                           SUNPSetupFn Psetup,
                                           SUNPSolveFn Psolve)
{
  /* set function pointers to integrator-supplied Psetup and PSolve
     routines and data, and return with success */
  BGMR_CONTENT(S)->Psetup = Psetup;
  BGMR_CONTENT(S)->Psolve = Psolve;
  BGMR_CONTENT(S)->PData  = PData;
  return SUN_SUCCESS;
}

int SUNLinSolSetup_BGMR(SUNLinearSolver S, [[maybe_unused]] SUNMatrix A)
{
  int status = SUN_SUCCESS;

  /* Set shortcuts to BGMR memory structures */
  SUNPSetupFn Psetup = BGMR_CONTENT(S)->Psetup;
  void* PData        = BGMR_CONTENT(S)->PData;

  /* no solver-specific setup is required, but if user-supplied
     Psetup routine exists, call that here */
  if (Psetup != NULL)
  {
    status = Psetup(PData);
    if (status != 0)
    {
      // options, SUNLS_PSET_FAIL_UNREC or SUNLS_PSET_FAIL_REC
      // NOTE: schur_psetup only returns 0 so this isnt a possible branch currently
      return (SUNLS_PSET_FAIL_UNREC);
    }
  }

  /* return with success */
  return SUN_SUCCESS;
}

int SUNLinSolSolve_BGMR(SUNLinearSolver S, [[maybe_unused]] SUNMatrix A,
                         N_Vector x, N_Vector b, sunrealtype delta)
{
  auto* content = BGMR_CONTENT(S);
  return content->solve(x, b, delta);
}

int SUNLinSolNumIters_BGMR(SUNLinearSolver S)
{
  return (BGMR_CONTENT(S)->numiters);
}

sunrealtype SUNLinSolResNorm_BGMR(SUNLinearSolver S)
{
  return (BGMR_CONTENT(S)->resnorm);
}

sunindextype SUNLinSolLastFlag_BGMR(SUNLinearSolver S)
{
  return (BGMR_CONTENT(S)->last_flag);
}

SUNErrCode SUNLinSolFree_BGMR(SUNLinearSolver S)
{ 
  if (S == nullptr) {
    return SUN_SUCCESS;
  }
  auto* content = BGMR_CONTENT(S);
  delete content;
  S->content = nullptr;
  SUNLinSolFreeEmpty(S);
  return SUN_SUCCESS;
}

} // namespace TChem