Running the LiLF pipeline
===============================

The 'standard' running protocol for LiLF contains five separate steps. Each step should be run one after another, and is designed such that the steps can easily be joined to form a singular pipeline.

LiLF is used to calibrate one or several LOFAR measurement sets of a specific target pointing. This target pointing can be observed over several nights.
Each observations should contain both the target pointing, as well as a calibrator observation (simultaneously). Currently, there is support for the following calibrators:

* 3C196
* 3C380
* 3C295 (not decameter/IS)
* 3C147 (only for international stations) 
* 3C 48 (only for international stations)


LOFAR_cal
^^^^^^^^^^^^^^^^^^^^^^^^^
LOFAR_cal takes 
