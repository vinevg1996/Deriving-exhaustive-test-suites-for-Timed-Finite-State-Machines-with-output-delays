# Deriving exhaustive test suites for Timed Finite State Machines with output delays

## Reproducing the Experiments

To reproduce the experiments with **random TFSMs** or **case-study TFSMs**, run:

```bash
python3 main.py
```

Then follow the instructions displayed by the program.

### Random TFSMs

Random machines with \(n = 5, 6, \ldots, 20\) states are located in:

```text
random_fsms/fsms/
```

### Case-Study FSMs and TFSMs

The FSMs and TFSMs used in the case studies are located in:

```text
case_studies/ansible/
case_studies/sdn/
```

Additional FSMs are also available in:

```text
random_fsms/fsms/
```

### Deriving the Ansible Specification TFSM

To derive the specification TFSM that models Ansible, run the playbooks used for time evaluation. These playbooks are located in:

```text
case_studies/ansible/playbooks_for_time_evaluation/
```
