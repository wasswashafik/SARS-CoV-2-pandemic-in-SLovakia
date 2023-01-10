# SK-COVID-19

 Data and MATLAB analysis for the SARS-CoV-2 pandemic in SLovakia.

I'm a control engineer, so take this with a pinch of salt. Data hand-entered based on the local Public Health Authority of the Slovak Republic (http://www.uvzsr.sk/en/) data. (Even though it would be their job to share the data online). The extrapolation is a simple exponential data fit with 95% confidence intervals at the moment. The crude estimation of the total cases is based on the study (https://annals.org/aim/fullarticle/2762808/incubation-period-coronavirus-disease-2019-covid-19-from-publicly-reported) claiming a mean infection-to-symptom shift of 5 days, where the confirmed case prediction is used. 

If you have an idea to improve this, let me know - or even better - Fork and open a PR.

# FAQ

- Can you predict when will it end and how many cases will we have? At the moment no. The closest we can have to such a prediction is when the infections slow down and we pass the exponential stage.
- Do your predictions include the aggessive prevention measures taken by the Slovak state? No. All I can project from is the data we have.
- Does the model account for the number of testing in relation to population size? No. I'm afraid it is not so easy to account for these factors.
