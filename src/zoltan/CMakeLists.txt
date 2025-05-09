project(zoltan1 C)

file(GLOB zoltan1_src src/all/all_allo.c
          src/coloring/coloring.c src/coloring/color_test.c
          src/coloring/bucket.c src/coloring/g2l_hash.c
          src/graph/graph.c src/ha/divide_machine.c
          src/ha/get_processor_name.c src/ha/ha_ovis.c
          src/hier/hier.c src/hier/hier_free_struct.c
          src/hsfc/hsfc_box_assign.c src/hsfc/hsfc.c
          src/hsfc/hsfc_hilbert.c
          src/hsfc/hsfc_point_assign.c src/lb/lb_balance.c
          src/lb/lb_box_assign.c src/lb/lb_copy.c
          src/lb/lb_eval.c src/lb/lb_free.c
          src/lb/lb_init.c src/lb/lb_invert.c
          src/lb/lb_migrate.c src/lb/lb_part2proc.c
          src/lb/lb_point_assign.c src/lb/lb_remap.c
          src/lb/lb_set_fn.c src/lb/lb_set_method.c
          src/lb/lb_set_part_sizes.c
          src/matrix/matrix_build.c
          src/matrix/matrix_distribute.c
          src/matrix/matrix_operations.c
          src/matrix/matrix_sym.c src/matrix/matrix_utils.c
          src/order/order.c src/order/order_struct.c
          src/order/order_tools.c src/order/hsfcOrder.c
          src/order/perm.c src/par/par_average.c
          src/par/par_bisect.c src/par/par_median.c
          src/par/par_median_randomized.c
          src/par/par_stats.c src/par/par_sync.c
          src/par/par_tflops_special.c
          src/params/assign_param_vals.c
          src/params/bind_param.c src/params/check_param.c
          src/params/free_params.c src/params/key_params.c
          src/params/print_params.c src/params/set_param.c
          src/tpls/build_graph.c src/tpls/postprocessing.c
          src/tpls/preprocessing.c src/tpls/scatter_graph.c
          src/tpls/third_library.c src/tpls/verify_graph.c
          src/tpls/parmetis_interface.c src/phg/phg_build.c
          src/phg/phg_build_calls.c src/phg/phg.c
          src/phg/phg_lookup.c src/phg/phg_verbose.c
          src/phg/phg_coarse.c src/phg/phg_comm.c
          src/phg/phg_distrib.c src/phg/phg_gather.c
          src/phg/phg_hypergraph.c src/phg/phg_match.c
          src/phg/phg_order.c src/phg/phg_parkway.c
          src/phg/phg_patoh.c src/phg/phg_plot.c
          src/phg/phg_rdivide.c src/phg/phg_refinement.c
          src/phg/phg_scale.c src/phg/phg_serialpartition.c
          src/phg/phg_util.c src/phg/phg_tree.c
          src/phg/phg_Vcycle.c src/rcb/box_assign.c
          src/rcb/create_proc_list.c src/rcb/inertial1d.c
          src/rcb/inertial2d.c src/rcb/inertial3d.c
          src/rcb/point_assign.c src/rcb/rcb_box.c
          src/rcb/rcb.c src/rcb/rcb_util.c
          src/rcb/rib.c src/rcb/rib_util.c
          src/rcb/shared.c src/reftree/reftree_build.c
          src/reftree/reftree_coarse_path.c
          src/reftree/reftree_hash.c
          src/reftree/reftree_part.c
          src/simple/block.c
          src/simple/cyclic.c src/simple/random.c
          src/timer/timer_params.c
          src/Utilities/Communication/comm_exchange_sizes.c
          src/Utilities/Communication/comm_invert_map.c
          src/Utilities/Communication/comm_do.c
          src/Utilities/Communication/comm_do_reverse.c
          src/Utilities/Communication/comm_info.c
          src/Utilities/Communication/comm_create.c
          src/Utilities/Communication/comm_resize.c
          src/Utilities/Communication/comm_sort_ints.c
          src/Utilities/Communication/comm_destroy.c
          src/Utilities/Communication/comm_invert_plan.c
          src/Utilities/Timer/zoltan_timer.c
          src/Utilities/Timer/timer.c
          src/Utilities/DDirectory/DD_Memory.c
          src/Utilities/DDirectory/DD_Find.c
          src/Utilities/DDirectory/DD_Destroy.c
          src/Utilities/DDirectory/DD_Set_Neighbor_Hash_Fn3.c
          src/Utilities/DDirectory/DD_Remove.c
          src/Utilities/DDirectory/DD_Create.c
          src/Utilities/DDirectory/DD_Update.c
          src/Utilities/DDirectory/DD_Stats.c
          src/Utilities/DDirectory/DD_Hash2.c
          src/Utilities/DDirectory/DD_Print.c
          src/Utilities/DDirectory/DD_Set_Neighbor_Hash_Fn2.c
          src/Utilities/DDirectory/DD_Set_Hash_Fn.c
          src/Utilities/DDirectory/DD_Set_Neighbor_Hash_Fn1.c
          src/Utilities/Memory/mem.c
          src/Utilities/shared/zoltan_align.c
          src/Utilities/shared/zoltan_id.c src/zz/zz_coord.c
          src/zz/zz_gen_files.c src/zz/zz_hash.c
          src/zz/murmur3.c src/zz/zz_map.c
          src/zz/zz_heap.c src/zz/zz_init.c
          src/zz/zz_obj_list.c src/zz/zz_rand.c
          src/zz/zz_set_fn.c src/zz/zz_sort.c
          src/zz/zz_struct.c src/zz/zz_back_trace.c
          src/zz/zz_util.c
)

file(GLOB zoltan_inc src/* src/Utilities/shared src/Utilities/Timer
     LIST_DIRECTORIES true)
add_library(zoltan1 ${zoltan1_src})
target_include_directories(zoltan1 PUBLIC ${zoltan_inc} ${MPI_C_INCLUDE_DIRS})

INSTALL(TARGETS zoltan1
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)
