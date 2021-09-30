THis folder use the BRCM toolbox to generate the one room and height rooms model. Also includes the reduction order of the height rooms model.

* Physical parameters can be modified in the excel file in the folder "Data_bench_coloc" (8 rooms model) and "Data_bench_solo" (1 room model) see BRCM documentation for details
* the matlab file bench_1_room.m and bench_8_rooms.m use the BRCM toolbox to generate the model.
* Model reduction is implementation in the dedicated folder for the 8 rooms model
* Some test of MPC basic MPC control are implemented in the loop files. This file use all the file "init..." and "plot..".

