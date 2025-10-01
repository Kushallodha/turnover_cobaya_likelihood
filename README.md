# Turnover Likelihood

A Cobaya likelihood for the DESI turnover observable DV(z) * k_TO(z).

## Structure

- Example Cobaya input: `turnover_likelihood/examples/cobaya_turnover_y1.yaml`

## Usage
Run the example:

```bash
cobaya-run turnover_likelihood/examples/cobaya_turnover_y1.yaml
```



## Credits
- Adapted from the turnover likelihood implementation in [desilike](https://github.com/cosmodesi/desilike) developed by Benedict Bahr-Kalus, Arnaud de Mattia and  David Parkinson. 

# Citations

If you are using the likelihood, please consider citing following works: 
```
@article{DESI:2025euz,
    author = "Bahr-Kalus, B. and others",
    collaboration = "DESI",
    title = "{Model-independent measurement of the matter-radiation equality scale in DESI 2024}",
    eprint = "2505.16153",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    reportNumber = "FERMILAB-PUB-25-0369-PPD",
    doi = "10.1103/yqm1-ybbv",
    journal = "Phys. Rev. D",
    volume = "112",
    number = "6",
    pages = "063553",
    year = "2025"
}

@article{Bahr-Kalus:2023ebd,
    author = "Bahr-Kalus, Benedict and Parkinson, David and Mueller, Eva-Maria",
    title = "{Measurement of the matter-radiation equality scale using the extended baryon oscillation spectroscopic survey quasar sample}",
    eprint = "2302.07484",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    doi = "10.1093/mnras/stad1867",
    journal = "Mon. Not. Roy. Astron. Soc.",
    volume = "524",
    number = "2",
    pages = "2463--2476",
    year = "2023",
    note = "[Erratum: Mon.Not.Roy.Astron.Soc. 526, 3248--3249 (2023)]"
}
``` 