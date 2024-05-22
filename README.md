# adzuki.bean
The codes used in the adzuki bean project

For most of the piplines, "qsub_subroutine.pl" file is required.

This file contains subroutines for different script files, including scripts for the PBS/PBS pro and Slurm server systems.

Please put "qsub_subroutine.pl" in your "$HOME" folder or "$HOME/softwares" folder.

The "create_job.pl" can independently generate scripts for job sending on server systems (qsub_subroutine.pl is also required).

Please modified IP address of the server in "qsub_subroutine.pl":

PBS server:
line 95

PBS pro server:
line 111

Slurm server:
line 98
