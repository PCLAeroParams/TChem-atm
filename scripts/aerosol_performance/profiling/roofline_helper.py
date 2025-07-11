import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#perf_scaling_factor = 1e11 # Nsight Compute default, 100 GFLOP/s
perf_scaling_factor = 1e9 # 1 GFLOP/s

def calculateRooflineStats(df, verbose=False):
    total_flops = 0
    total_dram_bytes = 0
    arith_intensity_arr = np.zeros(df.shape[0])
    performance_arr = np.zeros(df.shape[0])
    kernel_names = np.zeros(df.shape[0]).astype(str)
    timing_arr = np.zeros(df.shape[0])

    flop_count_units = {'smsp__sass_thread_inst_executed_op_dadd_pred_on.sum.per_cycle_elapsed': '', 
                     'smsp__sass_thread_inst_executed_op_dmul_pred_on.sum.per_cycle_elapsed': '',
                     'derived__smsp__sass_thread_inst_executed_op_dfma_pred_on_x2': ''}
    mem_units = {'dram__bytes.sum.per_second': ''}
    time_units = {'smsp__cycles_elapsed.avg.per_second': ''}

    for key in flop_count_units.keys():
        key_units = df[key][0]
        #print(key, key_units)
        flop_count_units[key] = key_units
    for key in mem_units.keys():
        key_units = df[key][0]
        #print(key, key_units)
        mem_units[key] = key_units
    for key in time_units.keys():
        key_units = df[key][0]
        #print(key, key_units)
        time_units[key] = key_units

    # Ensure FLOP units are either inst/cycle or inst
    if not all(units in ['inst/cycle', 'inst'] for units in flop_count_units.values()):
        raise AttributeError(f'Unexpected unit encountered in FLOP counts: {flop_count_units.values()}')

    # Determine units of memory bandwidth
    if mem_units['dram__bytes.sum.per_second'] == 'Gbyte/second':
        mem_conversion_factor = 1e9 # convert GB/s to B/s
    elif mem_units['dram__bytes.sum.per_second'] == 'Tbyte/second':
        mem_conversion_factor = 1e12 # convert TB/s to B/s
    else:
        raise AttributeError(f'Unexpected unit encountered in memory bandwith: {mem_units['dram__bytes.sum.per_second']}')

    # Determine clock speed units
    if time_units['smsp__cycles_elapsed.avg.per_second'] == 'cycle/nsecond':
        time_conversion_factor = 1e9 # convert cycle/nsecond (GHz) to cycle/second (Hz)
    else:
        raise AttributeError(f'Unexpected unit encountered in clockspeed: {time_units['smsp__cycles_elapsed.avg.per_second']}')

    for i in range(1, df.shape[0]):
        kernel = df.iloc[i, :]
        name = kernel['Kernel Name']

        mangle_chars = [
            'void Kokkos::',
            'Impl::',
            'cuda_parallel_launch_local_memory<Kokkos::Impl::ParallelFor<void sundials::kokkos::impl::',
            '<sundials::kokkos::Vector<Kokkos::Cuda, Kokkos::CudaSpace>>',
            '(double, _generic_N_Vector *)::',
            '[lambda(unsigned int) (instance 1)],',
            'Kokkos::RangePolicy<Kokkos::Cuda>, Kokkos::Cuda>>(T1)',
            'cuda_parallel_launch_local_memory<Kokkos::ParallelFor<void sundials::kokkos::impl::',
            '_generic_N_Vector *, _generic_N_Vector *)::', 
            '<unnamed>::',
            '(double, _generic_N_Vector *, double,', 
            'cuda_parallel_launch_local_memory<Kokkos::ParallelReduce<Kokkos::CombinedFunctorReducer<double sundials::kokkos::impl::',
            '([lambda(unsigned int, double &) (instance 1)], Kokkos::FunctorAnalysis<Kokkos::FunctorPatternInterface::REDUCE, Kokkos::RangePolicy<Kokkos::Cuda>, double sundials::kokkos::impl::N_VDotProd_Kokkos([lambda(unsigned int, double &) (instance 1)], double>::Reducer, void>, ', 
            '([lambda(unsigned int, double &) (instance 1)], Kokkos::FunctorAnalysis<Kokkos::FunctorPatternInterface::REDUCE, Kokkos::RangePolicy<Kokkos::Cuda>, double sundials::kokkos::impl::N_VWSqrSumLocal_Kokkos([lambda(unsigned int, double &) (instance 1)], double>::Reducer, void>, ',
            'cuda_parallel_launch_constant_memory<Kokkos::ParallelFor<TChem::', 
            ', Kokkos::TeamPolicy<Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static>>, Kokkos::Cuda>>()',
            '::[lambda(unsigned int, double &) (instance 1)], Kokkos::FunctorAnalysis<Kokkos::FunctorPatternInterface::REDUCE, Kokkos::RangePolicy<Kokkos::Cuda>, ',
            'Kokkos::Max<double, Kokkos::HostSpace>, double>::Reducer, void>, ',
            'Kokkos::Min<double, Kokkos::HostSpace>, double>::Reducer, void>, ',
            '(double,', '(_generic_N_Vector *)', '(_generic_N_Vector *,',
            '<Kokkos::ParallelFor<main:: Kokkos::RangePolicy<Kokkos::Cuda>, Kokkos::Cuda>>()', 'int *)', '(', ')', 
            ' '
        ]
        for string in mangle_chars:
            name = name.replace(string, '')
        
        kernel_names[i] = name
        achieved_work = (
            float(kernel['smsp__sass_thread_inst_executed_op_dadd_pred_on.sum.per_cycle_elapsed']) + 
            float(kernel['smsp__sass_thread_inst_executed_op_dmul_pred_on.sum.per_cycle_elapsed']) + 
            float(kernel["derived__smsp__sass_thread_inst_executed_op_dfma_pred_on_x2"])
        )*float(kernel['smsp__cycles_elapsed.avg.per_second'])*time_conversion_factor
        achieved_trafic = float(kernel["dram__bytes.sum.per_second"])*mem_conversion_factor

        total_flops += achieved_work
        total_dram_bytes += achieved_trafic 

        arithmetic_intensity = achieved_work / achieved_trafic
        performance = achieved_work

        #if i == 11 and verbose:
        if i == 20 and verbose:
            print(f"{i:.0f}: {name}")
            #print("--- Aggregation Complete ---")
            #print(f"Total Kernels (rows in CSV): {len(df_data)}")
            print(f"\tFP64 FLOPs: {achieved_work:,.0f}")
            print(f"\tDRAM Bytes Transferred: {achieved_trafic:,.0f} B")
            #print(f"Total Kernel Duration: {total_duration_s:.6f} s")
            #print("-" * 30)

            print(f"\tArithmetic Intensity: {arithmetic_intensity:,.2f} [FLOP/Byte]")
            print(f"\tPerformance: {performance:,.2f} [FLOP/s]")

        arith_intensity_arr[i] = arithmetic_intensity
        performance_arr[i] = performance/perf_scaling_factor
        timing_arr[i] = kernel['gpu__time_duration.avg']
    #total_arithmetic_intensity = total_flops / total_dram_bytes
    #total_performance = total_flops
    return arith_intensity_arr, performance_arr, kernel_names, timing_arr

def plotRoofline(arith_intensity_arr, performance_arr, kernel_names, timing_arr=np.array([]), **kwargs):
    
    if not kwargs.get('ax', None):
        fig, ax = plt.subplots(1, 1)
    else:
        ax = kwargs.get('ax')

    ylims = kwargs.get('ylim', (1e-3, 1e4))
    xlims = kwargs.get('xlim', (1e-1, 1e4))
    ax.set_ylim(ylims[0], ylims[1])
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(xlims[0], xlims[1])
    
    H100_fp64_peak_work = 25817357584788
    peak_traffic = 4022783332471.05

    arith_intes_at_peak_perf = H100_fp64_peak_work/peak_traffic
    arith_intens_arr = np.logspace(np.log10(xlims[0]), np.log10(arith_intes_at_peak_perf), 100)
    traffic_line = peak_traffic*arith_intens_arr/perf_scaling_factor

    ax.plot((arith_intes_at_peak_perf, xlims[1]), (H100_fp64_peak_work/perf_scaling_factor, 
                                              H100_fp64_peak_work/perf_scaling_factor), c='#172864')
    ax.plot(arith_intens_arr, traffic_line, c='#172864')

    #print(timing_arr.size)

    if timing_arr.size == 0:
        for i, name in enumerate(np.unique(kernel_names)[1:], 1):

            mask = kernel_names == name
            if 'N_' in name:
                c='k'
            elif 'RHS' in name:
                c='r'
            else:
                c='b'

            ax.scatter(arith_intensity_arr[mask], performance_arr[mask], c=c, s=10)

    if perf_scaling_factor == 1e11:
        perf_units = '100 GLOP/s'
    elif perf_scaling_factor == 1e9:
        perf_units = 'GLOP/s'
    else:
        raise AttributeError(f'Invalid performance scaling factor: {perf_scaling_factor}')
    fontsize = kwargs.get('fontsize', 10)
    ax.set_xlabel('Arithmetic Intensity (FLOP/byte)', fontsize=fontsize)
    ax.set_ylabel(f'Performance ({perf_units})', fontsize=fontsize)
    ax.tick_params(axis='both', which='major', labelsize=fontsize)

    if timing_arr.size != 0:
        total_arith_intens = ((arith_intensity_arr*timing_arr)/(timing_arr.sum())).sum() #arith_intens_arr.sum()/full_time
        total_performance = ((performance_arr*timing_arr)/(timing_arr.sum())).sum() #performance_arr.sum()/full_time
        
        kwargs_remove = ['ax', 'xlim', 'ylim', 'fontsize']
        for arg in kwargs_remove:
            if arg in kwargs:
                kwargs.pop(arg)

        ax.scatter(total_arith_intens, total_performance, **kwargs)
        #print(total_arith_intens, total_performance)
    #ax.axvline(6.42, ymin=0, ymax=1)

def formatPlotGrid(ax, **kwargs):
    """Helper function for formatting plot grid
    """
    major_lw = kwargs.get('major_lw', 1)
    minor_lw = kwargs.get('minor_lw', 1)
    tick_labelsize = kwargs.get('tick_labelsize', 10)
    ax.grid(which = "major", linewidth = major_lw, axis='y', ls="dotted", dashes=(major_lw,6), c='#414141', alpha=.5)
    ax.grid(which = "minor", linewidth = minor_lw, axis='y', ls="dotted", dashes=(minor_lw,6), c='white')
    ax.grid(which = "minor", linewidth = minor_lw, axis='x', ls="dotted", dashes=(minor_lw,6), c='#414141')
    ax.grid(which = "major", linewidth = major_lw, axis='x', ls="dotted", dashes=(major_lw,6), c='#414141')
    ax.tick_params(axis='both', labelsize=tick_labelsize, which='major', direction='in', top=True, right=True, bottom=True, left=True)
    ax.tick_params(axis='both', which='minor',direction='in',top=True, right=True, bottom=True, left=True)
