# Risk Factors Associated with Low Infant Birth Weight

This version of the data set has been deprecated and will not be
supported in future versions. Please use
[`Birthwt`](https://pbreheny.github.io/grpreg/reference/Birthwt.md)
instead.

## Usage

``` r
birthwt.grpreg
```

## Format

This data frame contains the following columns:

- `low` Indicator of birth weight less than 2.5kg

- `bwt` Birth weight in kilograms

- `age1,age2,age3` Orthogonal polynomials of first, second, and third
  degree representing mother's age in years

- `lwt1,lwt2,lwt3` Orthogonal polynomials of first, second, and third
  degree representing mother's weight in pounds at last menstrual period

- `white,black` Indicator functions for mother's race; "other" is
  reference group

- `smoke` smoking status during pregnancy

- `ptl1,ptl2m` Indicator functions for one or for two or more previous
  premature labors, respectively. No previous premature labors is the
  reference category.

- `ht` History of hypertension

- `ui` Presence of uterine irritability

- `ftv1,ftv2,ftv3m` Indicator functions for one, for two, or for three
  or more physician visits during the first trimester, respectively. No
  visits is the reference category.

## See also

[`Birthwt`](https://pbreheny.github.io/grpreg/reference/Birthwt.md)
