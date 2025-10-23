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


`LOFAR_preprocess <desc/preprocess>`
^^^^^^^^^^^^^^^^^^^^^^^^^
LOFAR_preprocess is the pipeline used to download data from the LOFAR `Long Term Archive (LTA) <https://lta.lofar.eu/>`_. Downloading data from the LTA requires an account with LOFAR/ASTRON. In order to download data from the LTA, the login credentials should be saved in a plain text file in your home directory (.wgetrc):
.. code-block::
  user=exampleUsername
  password=1234testPassword1234

Given that the login credentials are saved in plain text, it is STRONGLY encouraged not to use a password that you use for anything else.



LOFAR_cal
^^^^^^^^^^^^^^^^^^^^^^^^^
LOFAR_cal takes 
