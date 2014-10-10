OPTTYPE=bestfit
OPTTYPE=heating
VGLDISPLAY=:7.0

#JOBLIST = Study0035/0530 Study0023/0433 Study0023/0428 Study0023/0425 Study0030/0495 Study0030/0497 Study0030/0488 Study0030/0491 Study0030/0496 Study0030/0490 Study0017/0378 Study0018/0388 Study0018/0402 Study0018/0389 Study0018/0385 Study0029/0476 Study0029/0477 Study0025/0438 Study0025/0435 Study0025/0440 Study0025/0436 Study0028/0466 Study0028/0468 Study0028/0471 Study0052/0725 Study0052/0720 Study0026/0447 Study0026/0457 Study0026/0455 Study0026/0453 Study0026/0450 Study0026/0451 Study0057/0772 Study0057/0769 Study0022/0418 Study0022/0417 Study0021/0409 Study0021/0414 Study0021/0415 Study0054/0753 Study0054/0756 Study0053/0755 

JOBLIST = Study0023/0428 Study0030/0495 Study0030/0497 Study0030/0488 Study0030/0491 Study0030/0496 Study0030/0490 Study0017/0378 Study0018/0402 Study0018/0389 Study0018/0385 Study0029/0476 Study0029/0477 Study0025/0438 Study0025/0435 Study0025/0440 Study0025/0436 Study0028/0466 Study0028/0468 Study0028/0471 Study0026/0447 Study0026/0457 Study0026/0455 Study0026/0453 Study0026/0450 Study0026/0451 Study0022/0418 Study0022/0417 Study0021/0409 Study0021/0414 Study0021/0415


NUMGPU := "x x x x x"
LOWERCOUNT := "x"
UPPERCOUNT := $(NUMGPU)

%.gpu: 
	@echo "$@ $(words $(LOWERCOUNT)) $(words $(UPPERCOUNT)) "
	$(foreach var,$(wordlist  $(words $(LOWERCOUNT)), $(words $(UPPERCOUNT)),$(JOBLIST)),export GPUWORKDIR="optpp_pds/$(@:.gpu=)";echo $(var); dakota ./workdir/$(var)/opt/dakota_q_newton_$(OPTTYPE).in > ./workdir/$(var)/opt/dakota_q_newton_$(OPTTYPE).in.log 2>&1; export DISPLAY=$(VGLDISPLAY); python ./brainsearch.py --run_min ./workdir/$(var)/opt/optpp_pds.$(OPTTYPE) >> ./workdir/$(var)/opt/dakota_q_newton_$(OPTTYPE).in.log 2>&1;)
	#$(foreach var,$(wordlist  $(words $(LOWERCOUNT)), $(words $(UPPERCOUNT)),$(JOBLIST)),export GPUWORKDIR="optpp_pds/$(@:.gpu=)"; python ./brainsearch.py --run_min ./workdir/$(var)/opt/optpp_pds.$(OPTTYPE) >> ./workdir/$(var)/opt/dakota_q_newton_$(OPTTYPE).in.log 2>&1;)
#	@for job in $(wordlist  $(words $(LOWERCOUNT)), $(words $(UPPERCOUNT)),$(JOBLIST)) ;do\
#		echo "$$job gpu $(@:.gpu=) ./workdir/$$job/opt/dakota_q_newton_$(OPTTYPE).in.log ";\
#		export GPUWORKDIR="optpp_pds/$(@:.gpu=)";\
#		dakota ./workdir/$$job/opt/dakota_q_newton_$(OPTTYPE).in > ./workdir/$$job/opt/dakota_q_newton_$(OPTTYPE).in.log 2>&1;\
#		;\
#	done
	$(eval LOWERCOUNT += $(NUMGPU) )
	$(eval UPPERCOUNT += $(NUMGPU) )

all: 0.gpu 1.gpu 2.gpu 3.gpu 4.gpu 5.gpu 
