
### k12  (serial baseline T = 2.8926 s)

| variant | p | median time (s) | runs | speedup S | efficiency E |
|---|---|---|---|---|---|
| atomic_B | 1 | 2.9508 | 3 | 0.98 | 0.980 |
| atomic_B | 2 | 1.9673 | 3 | 1.47 | 0.735 |
| atomic_B | 4 | 1.6341 | 3 | 1.77 | 0.443 |
| private_C | 1 | 3.5842 | 3 | 0.81 | 0.807 |
| private_C | 2 | 2.7261 | 3 | 1.06 | 0.531 |
| private_C | 4 | 2.9018 | 3 | 1.00 | 0.249 |

### k8  (serial baseline T = 0.3390 s)

| variant | p | median time (s) | runs | speedup S | efficiency E |
|---|---|---|---|---|---|
| atomic_B | 1 | 0.7786 | 3 | 0.44 | 0.435 |
| atomic_B | 2 | 0.8221 | 3 | 0.41 | 0.206 |
| atomic_B | 4 | 0.7365 | 3 | 0.46 | 0.115 |
| private_C | 1 | 0.1633 | 3 | 2.08 | 2.077 |
| private_C | 2 | 0.0950 | 3 | 3.57 | 1.783 |
| private_C | 4 | 0.1249 | 3 | 2.71 | 0.678 |
