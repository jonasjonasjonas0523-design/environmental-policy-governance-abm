# v1.2 extended verification (v1.1 battery + V13-V17)

| test_id      | description                                                                                            | expected              | actual                                                                         | status   |
|:-------------|:-------------------------------------------------------------------------------------------------------|:----------------------|:-------------------------------------------------------------------------------|:---------|
| V1.1-battery | all 24 v1.1 items on the v1.2 build                                                                    | 24 PASS               | tau* = 0.093267  cap correction = 26.8%  FAILURES = 0                          | PASS     |
| V13          | Ann(x,0)=0; Ann(.9025,8); Ann(.9025,9)                                                                 | 0; 5.742291; 6.182417 | 0; 5.742291; 6.182417                                                          | PASS     |
| V14          | t=T: Ann=0, participate iff p*nu_i>1                                                                   | 52 participants       | 52                                                                             | PASS     |
| V15          | M nondecreasing in h and remaining horizon (m fixed)                                                   | monotone              | monotone                                                                       | PASS     |
| V16          | byte-equal frozen tau* everywhere; seed-424242 re-solve exact; indep quad <=5e-4; retired check absent | 0.09326697147369385   | arch==True, run==True, resolve==True, |ind-diff|=1.25e-04, retired absent=True | PASS     |
| V17          | ABM voluntary decisions == frozen closed-form inequality on probe grid                                 | 80 tick-masks equal   | 80 checked, all equal = True                                                   | PASS     |
