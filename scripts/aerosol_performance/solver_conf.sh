rhss_params(){
  export numerical_params="--do_rhs=true \
                           --do_jac=true"
}

rhs_params(){
  export numerical_params="--do_rhs=true \
                           --do_jac=false"
}

jac_params(){
  export numerical_params="--do_rhs=false \
                           --do_jac=true"
}

trbdf_params(){
  t_iterPerInt=1
  min_dt='1e-3'
  max_dt='1'
  atol_t='1e-12'
  tol_time='1e-3'
  export numerical_params="--tol-time=$tol_time \
                    --time-iterations-per-interval=$t_iterPerInt \
                    --dtmin=$min_dt \
                    --jacobian-interval=4 \
		    --dtmax=$max_dt \
                    --atol-time=${atol_t} \
                    --atol-newton=1e-16 \
                    --rtol-newton=1e-8" 

}

sundials_gmres_params(){
min_dt='1'
atol_t='1e-16'
rtol_time='1e-8'
export numerical_params="--rtol-time=${rtol_time} \
                    --dtmin=${min_dt} \
                    --solver_type=1 \
		    --atol-time=${atol_t} "
}

sundials_dense_params(){
min_dt='1'
atol_t='1e-16'
rtol_time='1e-8'
export numerical_params="--rtol-time=${rtol_time} \
                    --dtmin=${min_dt} \
                    --solver_type=0 \
                    --atol-time=${atol_t} "
}

