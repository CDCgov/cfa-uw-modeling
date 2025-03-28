### NSFG 2017-2019
### Import, Clean

# nolint start

library(here)
library(tidyverse)

# Female Dataset -----------#################
dat <- readRDS(here("data", "surveys", "2017_2019_FemRespData.rds"))
colnames(dat) <- tolower(colnames(dat))
fem_vars_wide <- readxl::read_xlsx(here(
  "data", "surveys",
  "nsfg_vars_females_2017_2019.xlsx"
))
fvars <- tolower(fem_vars_wide$Variable)
fvars <- c(fvars, "secu", "sest") # add clustering vars
d <- dat |> select(all_of(fvars))

## Define Relationships for Current Partners & Summary Degree Vars ------
d <- d |>
  mutate(
    curr1 = ifelse(pcurrnt == 1, 1, ifelse(pcurrnt == 5, 0, NA)),
    curr2 = ifelse(pcurrnt2 == 1, 1, ifelse(pcurrnt2 == 5, 0, NA)),
    curr3 = ifelse(pcurrnt3 == 1, 1, ifelse(pcurrnt3 == 5, 0, NA)),
    hadsex = ifelse(hadsex == 2, 0, hadsex),
    sexonce = ifelse(sexonce == 2, 0, sexonce),
    cohever = ifelse(cohever == 2, 0, cohever),
    evmarcoh = ifelse(evmarcoh == 2, 0, evmarcoh),
    condvag = ifelse(condvag == 1, 1, ifelse(condvag == 5, 0, NA)),
    condoral = ifelse(condfell == 1, 1, ifelse(condfell == 5, 0, NA)),
    condanal = ifelse(condanal == 1, 1, ifelse(condanal == 5, 0, NA)),
    condsex = ifelse(condsexl == 1, 1, ifelse(condsexl == 5, 0, NA)),
    chlamtest = ifelse(chlamtst == 1, 1, ifelse(chlamtst == 5, 0, NA)),
    stdtest = ifelse(stdothr12 == 1 | chlamtest == 1, 1,
      ifelse(stdothr12 == 5 & chlamtest == 0, 0, NA)
    ),
    gon = ifelse(gon == 1, 1, ifelse(gon == 5, 0, NA)),
    chlam = ifelse(chlam == 1, 1, ifelse(chlam == 5, 0, NA)),
    herpes = ifelse(herpes == 1, 1, ifelse(herpes == 5, 0, NA)),
    hpv = ifelse(genwarts == 1, 1, ifelse(genwarts == 5, 0, NA)),
    syphilis = ifelse(syphilis == 1, 1, ifelse(syphilis == 5, 0, NA)),
    stdtrt = ifelse(stdtrt12 == 1, 1, ifelse(stdtrt12 == 5, 0, NA)),
    hivtst = ifelse(evhivtst == 1, 1, ifelse(evhivtst == 5, 0, NA)),
    infever = ifelse(infever == 1, 1, ifelse(infever == 5, 0, NA)),
    pidtreat = ifelse(pidtreat == 1, 1, ifelse(pidtreat == 5, 0, NA)),
    samesexany = ifelse(samesexany == 1, 1,
      ifelse(samesexany == 5, 0, NA)
    ),
    evrinject = ifelse(evrinject == 1, 1, ifelse(evrinject == 5, 0, NA)),
    malsht12 = ifelse(malsht12 == 1, 1, ifelse(malsht12 == 5, 0, NA)),
    prostfrq = ifelse(prostfrq == 1, 1, ifelse(prostfrq == 5, 0, NA)),
    dateapp = ifelse(dateapp == 1, 1, ifelse(dateapp == 5, 0, NA)),
    bisexprt = ifelse(bisexprt == 1, 1, ifelse(bisexprt == 5, 0, NA)),
    cond1 = ifelse(pswkcond1 == 1, 1, ifelse(pswkcond1 == 5, 0, NA)),
    once1 = ifelse(partdur1 == 997, 1, 0),
    once2 = ifelse(partdur2 == 997, 1, 0),
    once3 = ifelse(partdur3 == 997, 1, 0),
    fuse1 = ifelse(usefstp == 1, 1, ifelse(usefstp == 5, 0, NA)),
    fuse2 = ifelse(usefstp2 == 1, 1, ifelse(usefstp2 == 5, 0, NA)),
    fuse3 = ifelse(usefstp3 == 1, 1, ifelse(usefstp3 == 5, 0, NA)),
    luse1 = ifelse(uselstp == 1, 1, ifelse(uselstp == 5, 0, NA)),
    luse2 = ifelse(uselstp2 == 1, 1, ifelse(uselstp2 == 5, 0, NA)),
    luse3 = ifelse(uselstp3 == 1, 1, ifelse(uselstp3 == 5, 0, NA)),
    firstuse = ifelse(usefrsts == 1, 1, ifelse(usefrsts == 5, 0, NA)),
    hxagemar = ifelse(hxagemar > 97, NA, hxagemar),
    hxagemar2 = ifelse(hxagemar2 > 97, NA, hxagemar2),
    hxagemar3 = ifelse(hxagemar3 > 97, NA, hxagemar3),
    hisagecx = ifelse(hisagecx > 97, NA, hisagecx),
    cphisage = ifelse(cphisage > 97, NA, cphisage),
    p1yhsage = ifelse(p1yhsage > 97, NA, p1yhsage),
    p1yhsage2 = ifelse(p1yhsage2 > 97, NA, p1yhsage2),
    p1yhsage3 = ifelse(p1yhsage3 > 97, NA, p1yhsage3)
  )

# continuous vars
d <- d |>
  mutate(
    partdur1 = ifelse(partdur1 < 997, partdur1, NA),
    partdur2 = ifelse(partdur2 < 997, partdur2, NA),
    partdur3 = ifelse(partdur3 < 997, partdur3, NA),
    sex4wk = ifelse(pst4wksx < 997, pst4wksx, NA),
    cond4wk = ifelse(pswkcond2 < 997, pswkcond2, NA)
  )

# make change partdur 0 to partdur = 0_5
# duration can't be 0, but indicates <1 month
d <- d |>
  mutate(
    partdur1 = ifelse(partdur1 == 0, 0.5, partdur1),
    partdur2 = ifelse(partdur2 == 0, 0.5, partdur2),
    partdur3 = ifelse(partdur3 == 0, 0.5, partdur3)
  )

# relationship status of recent partners
# first pass - use p1yrelp (relationship to partner)
d <- d |>
  mutate(
    # r had more than one partner in last 12 months, get from p1yrelp
    # current or former husband
    rel1 = ifelse(p1yrelp %in% 1:7, 1,
      # current or former cohab
      ifelse(p1yrelp %in% 8:12, 2,
        ifelse(p1yrelp %in% c(13, 20), 3,
          # first sex partner who is not a husband or "other"
          NA
        )
      )
    ),
    rel2 = ifelse(p1yrelp2 %in% 1:7, 1,
      ifelse(p1yrelp2 %in% 8:12, 2,
        ifelse(p1yrelp2 %in% c(13, 20), 3, NA)
      )
    ),
    rel3 = ifelse(p1yrelp3 %in% 1:7, 1,
      ifelse(p1yrelp3 %in% 8:12, 2,
        ifelse(p1yrelp3 %in% c(13, 20), 3, NA)
      )
    )
  )

# if respondent is married to or cohabitating with first sexual partner, add to rel1 and mark as current
# # there are a few people (7) who are separated but they recently had sex with their spouse, I count this as
# a current sexual partner (NSFG computed var currprtt does not include these partners as current)
d <- d |>
  mutate(
    rel1 = ifelse(is.na(rel1) & whofstpr == 7, 1, rel1), # current marriage
    rel1 = ifelse(is.na(rel1) & whofstpr == 8, 2, rel1), # current cohab
    curr1 = ifelse(rel1 %in% 1:2 & whofstpr %in% 7:8 & parts1yr == 1, 1, curr1) # flag as current
  )

# if respondent is married or cohabitating (informally) but doesn't report sex in last year, still mark as current
d <- d |>
  mutate(curr1 = ifelse(is.na(curr1) & rel1 %in% 1:2 & rmarital %in% 1:2, 1, curr1))

# if rel still undefined, let's pull from relapt1-3 recodes (relationship at first sex)
d <- d |>
  mutate(
    rel1 = ifelse(is.na(rel1) & relatp1 %in% 1:2, 1, # married or engaged
      ifelse(is.na(rel1) & relatp1 == 3, 2, # cohabitating
        ifelse(is.na(rel1) & relatp1 %in% 4:8, 3, rel1)
      )
    ), # anything else
    rel2 = ifelse(is.na(rel2) & relatp2 %in% 1:2, 1,
      ifelse(is.na(rel2) & relatp2 == 3, 2,
        ifelse(is.na(rel2) & relatp2 %in% 4:8, 3, rel2)
      )
    ),
    rel3 = ifelse(is.na(rel3) & relatp3 %in% 1:2, 1,
      ifelse(is.na(rel3) & relatp3 == 3, 2,
        ifelse(is.na(rel3) & relatp3 %in% 4:8, 3, rel3)
      )
    )
  )

# finally, there are a few people people who have an undefined current partner who they haven't had sex with in last year
# mark them as rel1 == 3

d <- d |> mutate(
  rel1 = ifelse(is.na(rel1) & parts1yr == 0 & curr1 == 1, 3, rel1),
  rel2 = ifelse(is.na(rel2) & parts1yr == 0 & curr2 == 1, 3, rel2),
  rel3 = ifelse(is.na(rel3) & parts1yr == 0 & curr3 == 1, 3, rel3)
)

# make binary vars for type of rel
d <- d |> mutate(
  mar1 = ifelse(rel1 == 1, 1, 0),
  cohab1 = ifelse(rel1 == 2, 1, 0),
  casual1 = ifelse(rel1 == 3, 1, 0),
  mar2 = ifelse(rel2 == 1, 1, 0),
  cohab2 = ifelse(rel2 == 2, 1, 0),
  casual2 = ifelse(rel2 == 3, 1, 0),
  mar3 = ifelse(rel3 == 1, 1, 0),
  cohab3 = ifelse(rel3 == 2, 1, 0),
  casual3 = ifelse(rel3 == 3, 1, 0)
)


# indicator variables for active AND rel type
d <- d |> mutate(
  mar1_active = ifelse(mar1 == 1 & curr1 == 1, 1, 0),
  cohab1_active = ifelse(cohab1 == 1 & curr1 == 1, 1, 0),
  casual1_active = ifelse(casual1 == 1 & curr1 == 1, 1, 0),
  mar2_active = ifelse(mar2 == 1 & curr2 == 1, 1, 0),
  cohab2_active = ifelse(cohab2 == 1 & curr2 == 1, 1, 0),
  casual2_active = ifelse(casual2 == 1 & curr2 == 1, 1, 0),
  mar3_active = ifelse(mar3 == 1 & curr3 == 1, 1, 0),
  cohab3_active = ifelse(cohab3 == 1 & curr3 == 1, 1, 0),
  casual3_active = ifelse(casual3 == 1 & curr3 == 1, 1, 0)
)

# classify one-times as once == 1 and curr == 0 (had sex with once and does not consider to be a current partner)
d <- d |> mutate(
  inst1 = ifelse(once1 == 1 & curr1 == 0, 1, 0),
  inst2 = ifelse(once2 == 1 & curr2 == 0, 1, 0),
  inst3 = ifelse(once3 == 1 & curr3 == 0, 1, 0)
)

# summary current degree
d <- d |>
  rowwise() |>
  mutate(
    deg = sum(curr1, curr2, curr3, na.rm = T),
    deg_mar = sum(mar1_active, mar2_active, mar3_active, na.rm = T),
    deg_cohab = sum(cohab1_active, cohab2_active, cohab3_active, na.rm = T),
    deg_casual = sum(casual1_active, casual2_active, casual3_active, na.rm = T),
    deg_inst = sum(inst1, inst2, inst3, na.rm = T),
    sum_parts = sum(!is.na(curr1), !is.na(curr2), !is.na(curr3))
  )

d <- d |>
  mutate(deg_main = ifelse((deg_mar | deg_cohab) == 1, 1, 0))

d <- d |> mutate(partsdiff = parts1yr - 3) # max # of detailed partners

# assume all partners in last year beyond 3 were one-time
# (applies to very few participants)
d <- d |>
  mutate(deg_inst_high = ifelse(partsdiff > 0, deg_inst + partsdiff, deg_inst)) |>
  mutate(deg_inst_high = ifelse(is.na(deg_inst_high), 0, deg_inst_high))

## Define Age for current partnerships ------ ###################

# Respondent Age
d <- d |> mutate(age = ifelse(ager < 97, ager, NA))


# We have some information about partner and respondent age at the start of relationships
# Respondent Age
d <- d |> mutate(
  initage1 = ifelse(p1yrage < 97, p1yrage, NA), # age at init sex w/ recent partners (but not first ever partner)
  initage2 = ifelse(p1yrage2 < 97, p1yrage2, NA),
  initage3 = ifelse(p1yrage3 < 97, p1yrage3, NA),
  initage1 = ifelse(is.na(initage1) & p1yrelp == 13, vry1stag, initage1), # age at first sexual partnership
  initage2 = ifelse(is.na(initage2) & p1yrelp2 == 13, vry1stag, initage2),
  initage3 = ifelse(is.na(initage3) & p1yrelp3 == 13, vry1stag, initage3)
)

# Fill in some missing info based on current/former mar/cohs for special cases for most recent partner
d <- d |> mutate(
  initage1 = ifelse(is.na(initage1) & whofstpr == 7 & lifprtnr == 1, fmar1age, initage1), # current husband is first & only sexual partner
  initage1 = ifelse(is.na(initage1) & whofstpr == 8 & lifprtnr == 1, cpherage, initage1), # current cohab is first & only sexual partner
  initage1 = ifelse(is.na(initage1) & p1yrelp == 7 & prevhusb == 0, fmar1age, initage1), # current and only husband is most recent partner
  initage1 = ifelse(is.na(initage1) & p1yrelp == 8, cpherage, initage1), # current cohab
  initage1 = ifelse(is.na(initage1) & p1yrelp == 9, heragecx, initage1), # most recent sex is first former cohabitating partner
  initage1 = ifelse(is.na(initage1) & !is.na(rel1) & curr1 == 1, age - (partdur1 / 12), initage1), # For everyone else, get age at init from current age - duration
)

# R age for second and third partners (sometimes age at init for R and P unknown or refused)
d <- d |>
  mutate(
    initage2 = ifelse(is.na(initage2) & !is.na(rel1) & curr2 == 1, age - (partdur2 / 12), initage2),
    initage3 = ifelse(is.na(initage3) & !is.na(rel3) & curr3 == 1, age - (partdur3 / 12), initage3),
  )

# Partner Age at beginning of sexual relationships
# WHAT IS MISSING: age of partner if partner was R's first sexual partner and they are not married/cohab
# AND age of partner if R is under 18
# plus occaional other missingness, some unknown ages for casual partners

d <- d |>
  mutate(
    pinitage1 = p1yhsage,
    pinitage2 = p1yhsage2,
    pinitage3 = p1yhsage3
  )

# Fill in some missing info based on current/former mar/cohs
d <- d |>
  mutate(
    pinitage1 = ifelse(is.na(pinitage1) & rel1 == 1 &
      prevhusb == 0 & !is.na(hxagemar),
    hxagemar, pinitage1
    ), # current & first husband
    pinitage1 = ifelse(is.na(pinitage1) & p1yrelp == 8 & !is.na(cphisage),
      cphisage, pinitage1
    ), # current cohab
    pinitage1 = ifelse(is.na(pinitage1) & rel1 == 1 & prevhusb == 1 & !is.na(hxagemar2),
      hxagemar2, pinitage1
    ), # second husband
    pinitage1 = ifelse(is.na(pinitage1) & rel1 == 1 & prevhusb == 2 & !is.na(hxagemar3),
      hxagemar3, pinitage1
    ), # third husband
    pinitage1 = ifelse(is.na(pinitage1) & p1yrelp == 9 & !is.na(hisagecx),
      hisagecx, pinitage1
    ) # former cohab
  )


# Calc age differences and current partner age
d <- d |>
  mutate(
    agediff1 = abs(pinitage1 - initage1),
    page1 = round(pinitage1 + partdur1 / 12),
    agediff2 = abs(pinitage2 - initage2),
    page2 = round(pinitage2 + partdur2 / 12),
    agediff3 = abs(pinitage3 - initage3),
    page3 = round(pinitage3 + partdur3 / 12)
  )

## Define Race / Ethnicity ------
# surprisingly complete
# Respondent Race (Hispanic, NH White, NH Black, Other/Multiple)
d <- d |> mutate(race = hisprace2)

# Partner Race (for current partners)
d <- d |>
  mutate(
    prace1 = ifelse(fmarital == 1 & rel1 == 1 & prevhusb == 0, hsbnrace1, NA), # current and first husband
    prace1 = ifelse(fmarital == 1 & rel1 == 1 & prevhusb == 1, hsbnrace2, prace1), # current hub is 2nd ever hub
    prace1 = ifelse(fmarital == 4 & rel1 == 1 & prevhusb == 1, hsbnrace1, prace1), # divorced/separated husband
    prace1 = ifelse(fmarital == 1 & rel1 == 1 & prevhusb == 2, hsbnrace3, prace1), # 3rd hub
    prace1 = ifelse(fmarital == 4 & rel1 == 1 & prevhusb == 2, hsbnrace2, prace1),
    prace1 = ifelse(fmarital == 1 & rel1 == 1 & prevhusb == 3, hsbnrace4, prace1), # 4th hub
    prace1 = ifelse(fmarital == 4 & rel1 == 1 & prevhusb == 3, hsbnrace3, prace1),
    prace1 = ifelse(fmarital == 4 & rel1 == 1 & prevhusb == 4, hsbnrace4, prace1),
    prace1 = ifelse(is.na(prace1) & rmarital == 2 & rel1 == 2, curcohnrace, prace1), # current cohabitating partner
    prace1 = ifelse(is.na(prace1) & p1yrelp == 13 & curr1 == 1, fsexnrace, prace1), # race of first sexual partner who is still current
    prace1 = ifelse(is.na(prace1), p1ynrace1, prace1),
    prace2 = p1ynrace2,
    prace2 = ifelse(is.na(prace2) & p1yrelp2 == 13 & curr2 == 1, fsexnrace, prace2),
    prace2 = ifelse(is.na(prace2) & p1yrelp2 == 7 & fmarital == 1 & prevhusb == 0, hsbnrace1, prace2),
    prace2 = ifelse(is.na(prace2) & p1yrelp2 == 7 & fmarital == 1 & prevhusb == 1, hsbnrace2, prace2),
    prace2 = ifelse(is.na(prace2) & p1yrelp2 == 7 & fmarital == 4 & prevhusb == 1, hsbnrace1, prace2),
    prace2 = ifelse(is.na(prace2) & p1yrelp2 == 7 & fmarital == 4 & prevhusb == 3, hsbnrace3, prace2),
    prace2 = ifelse(is.na(prace2) & p1yrelp2 == 8, curcohnrace, prace2),
    prace3 = p1ynrace3,
    prace3 = ifelse(is.na(prace3) & p1yrelp3 == 13 & curr3 == 1, fsexnrace, prace3)
  )

## Condom Use ------- ########

d <- d |> mutate(
  # Overall in last month
  p_cond_month = NA,
  p_cond_month = ifelse(sex4wk > 1 & !is.na(cond4wk), cond4wk / sex4wk, p_cond_month), # sex >1, cond4wk not missing
  p_cond_month = ifelse(sex4wk > 1 & is.na(cond4wk), 0, p_cond_month), # sex >1, assume no condoms
  p_cond_month = ifelse(sex4wk == 1 & cond1 == 1, 1, p_cond_month), # one sex in last month, yes condom
  p_cond_month = ifelse(sex4wk == 1 & cond1 == 0, 0, p_cond_month), # one sex in last month, no condom
  p_cond_month = ifelse(p_cond_month >= 20, NA, p_cond_month), # one person had 5 sex and 100 condom? ignore
  # Condom use at first and last sex with recent partners
  firstevercond = NA,
  firstevercond = ifelse(firstuse == 1 & (mthfrsts1 == 4 | mthfrsts2 == 4), 1, 0),
  fsex1 = NA,
  fsex1 = ifelse(fuse1 == 1 & (fstmthp1 == 4 | fstmthp2 == 4), 1, 0),
  fsex2 = NA,
  fsex2 = ifelse(fuse2 == 1 & (fstmthp5 == 4 | fstmthp6 == 4), 1, 0),
  fsex3 = NA,
  fsex3 = ifelse(fuse3 == 1 & (fstmthp9 == 4 | fstmthp10 == 4), 1, 0),
  fsex1 = ifelse(is.na(fsex1) & p1yrelp == 13, firstevercond, fsex1),
  fsex2 = ifelse(is.na(fsex2) & p1yrelp2 == 13, firstevercond, fsex2),
  fsex3 = ifelse(is.na(fsex3) & p1yrelp3 == 13, firstevercond, fsex3),
  lsex1 = NA,
  lsex1 = ifelse(luse1 == 1 & (lstmthp1 == 4 | lstmthp2 == 4), 1, 0),
  lsex2 = NA,
  lsex2 = ifelse(luse2 == 1 & (lstmthp5 == 4 | lstmthp6 == 4), 1, 0),
  lsex3 = NA,
  lsex3 = ifelse(luse3 == 1 & (lstmthp9 == 4 | lstmthp10 == 4), 1, 0),
)

# is partner active within last month (or so, since condom question is "last 4 weeks")
# mutate(p1month = ifelse(cmintvw-cmlsx1 < 2, 1, 0),
#        p2month = ifelse(cmintvw-cmlsx2 < 2, 1, 0),
#        p3month = ifelse(cmintvw-cmlsx3 < 2, 1, 0)
# )


## Age Group (& 50 yr old adj) ------- #############
# make all 50-year olds 49.9 (they were 49 at screener interview,)
d <- d |>
  mutate(age = ifelse(age == 50, 49.9, age))

# make age groups
d$age_group <- cut(round(d$age), 7)

d <- d |> ungroup()

## Save out ------ ####

saveRDS(d, here("data", "surveys", "femdat.rds"))


# Male ----- #####
mdat <- readRDS(here("data", "surveys", "2017_2019_MaleData.rds"))
colnames(mdat) <- tolower(colnames(mdat))
male_vars_wide <- readxl::read_xlsx(here("data", "surveys", "nsfg_vars_males_2017_2019.xlsx"))
mvars <- tolower(male_vars_wide$Variable)
mvars <- c(mvars, "secu", "sest") # add clustering vars
m <- mdat |> select(all_of(mvars))

## Define Relationships for Current Partners & Summary Degree Vars ------ ###################
# binary vars
m <- m |>
  mutate(
    curr1 = ifelse(pxcurrprt == 1, 1, ifelse(pxcurrprt == 0, 0, NA)),
    curr2 = ifelse(pxcurrprt2 == 1, 1, ifelse(pxcurrprt2 == 0, 0, NA)),
    curr3 = ifelse(pxcurrprt3 == 1, 1, ifelse(pxcurrprt3 == 0, 0, NA)),
    hadsex = ifelse(hadsex == 2, 0, hadsex),
    sexonce = ifelse(sexonce == 2, 0, sexonce),
    cohever = ifelse(cohever == 2, 0, cohever),
    evmarcoh = ifelse(evmarcoh == 2, 0, evmarcoh),
    condvag = ifelse(condvag == 1, 1, ifelse(condvag == 5, 0, NA)),
    condoral = ifelse(condfell == 1, 1, ifelse(condfell == 5, 0, NA)),
    condanal = ifelse(condanal == 1, 1, ifelse(condanal == 5, 0, NA)),
    condsex = ifelse(condsexl == 1, 1, ifelse(condsexl == 5, 0, NA)),
    cond_sexonce = ifelse(p12mocono == 1, 1, ifelse(p12mocono == 5, 0, NA)),
    stdtest = ifelse(stdtst12 == 1, 1, ifelse(stdtst12 == 5, 0, NA)),
    gon = ifelse(gon == 1, 1, ifelse(gon == 5, 0, NA)),
    chlam = ifelse(chlam == 1, 1, ifelse(chlam == 5, 0, NA)),
    herpes = ifelse(herpes == 1, 1, ifelse(herpes == 5, 0, NA)),
    hpv = ifelse(genwarts == 1, 1, ifelse(genwarts == 5, 0, NA)),
    syphilis = ifelse(syphilis == 1, 1, ifelse(syphilis == 5, 0, NA)),
    stdtrt = ifelse(stdtrt12 == 1, 1, ifelse(stdtrt12 == 5, 0, NA)),
    hivtst = ifelse(evhivtst == 1, 1, ifelse(evhivtst == 5, 0, NA)),
    infever = ifelse(infever == 1, 1, ifelse(infever == 5, 0, NA)),
    samesexany = ifelse(samesexany == 1, 1, ifelse(samesexany == 5, 0, NA)),
    evrinject = ifelse(evrinject == 1, 1, ifelse(evrinject == 5, 0, NA)),
    femsht12 = ifelse(femsht12 == 1, 1, ifelse(femsht12 == 5, 0, NA)),
    prostfrq = ifelse(prostfrq == 1, 1, ifelse(prostfrq == 5, 0, NA)),
    once1 = ifelse(partdur1 == 997, 1, 0),
    once2 = ifelse(partdur2 == 997, 1, 0),
    once3 = ifelse(partdur3 == 997, 1, 0),
    firstusecwp = ifelse(cwpfuse == 1, 1, ifelse(cwpfuse == 5, 0, NA))
  )

# continuous vars
m <- m |>
  mutate(
    partdur1 = ifelse(partdur1 < 997, partdur1, NA),
    partdur2 = ifelse(partdur2 < 997, partdur2, NA),
    partdur3 = ifelse(partdur3 < 997, partdur3, NA),
    sex4wk = ifelse(sexfreq < 997, sexfreq, NA),
    cond4wk = ifelse(confreq < 997, confreq, NA),
    # cmlsx1 = ifelse(cmlsex < 9998, cmlsex, NA),
    # cmlsx2 = ifelse(cmlsex2 < 9998, cmlsex2, NA),
    # cmlsx3 = ifelse(cmlsex3 < 9998, cmlsex3, NA),
  )

# make change partdur 0 to partdur = 0_5
# duration can't be 0, but indicates <1 month
m <- m |>
  mutate(
    partdur1 = ifelse(partdur1 == 0, 0.5, partdur1),
    partdur2 = ifelse(partdur2 == 0, 0.5, partdur2),
    partdur3 = ifelse(partdur3 == 0, 0.5, partdur3)
  )

# relationship status of recent partners
# reclaim active status if partner is current wife or current cohab
m <- m |> mutate(
  curr1 = ifelse(p1relation %in% c(1, 2, 3), 1, curr1),
  curr2 = ifelse(p2relation %in% c(1, 2, 3), 1, curr2),
  curr3 = ifelse(p3relation %in% c(1, 2, 3), 1, curr3)
)

# relationship status at most recent sex with partner
m <- m |> mutate(
  rel1 = ifelse(p1relation %in% c(1, 2, 4), 1, # current or former wife
    ifelse(p1relation %in% c(3, 5), 2, # current or former cohab
      ifelse(p1relation %in% c(6, 8), 3, NA)
    )
  ), # "other"
  rel2 = ifelse(p2relation %in% c(1, 2, 4), 1,
    ifelse(p2relation %in% c(3, 5), 2,
      ifelse(p2relation %in% c(6, 8), 3, NA)
    )
  ),
  rel3 = ifelse(p3relation %in% c(1, 2, 4), 1,
    ifelse(p3relation %in% c(3, 5), 2,
      ifelse(p3relation %in% c(6, 8), 3, NA)
    )
  ),
  mar1 = ifelse(rel1 == 1, 1, 0),
  cohab1 = ifelse(rel1 == 2, 1, 0),
  casual1 = ifelse(rel1 == 3, 1, 0),
  mar2 = ifelse(rel2 == 1, 1, 0),
  cohab2 = ifelse(rel2 == 2, 1, 0),
  casual2 = ifelse(rel2 == 3, 1, 0),
  mar3 = ifelse(rel3 == 1, 1, 0),
  cohab3 = ifelse(rel3 == 2, 1, 0),
  casual3 = ifelse(rel3 == 3, 1, 0)
)

# indicator variables for active AND rel type
m <- m |> mutate(
  mar1_active = ifelse(mar1 == 1 & curr1 == 1, 1, 0),
  cohab1_active = ifelse(cohab1 == 1 & curr1 == 1, 1, 0),
  casual1_active = ifelse(casual1 == 1 & curr1 == 1, 1, 0),
  mar2_active = ifelse(mar2 == 1 & curr2 == 1, 1, 0),
  cohab2_active = ifelse(cohab2 == 1 & curr2 == 1, 1, 0),
  casual2_active = ifelse(casual2 == 1 & curr2 == 1, 1, 0),
  mar3_active = ifelse(mar3 == 1 & curr3 == 1, 1, 0),
  cohab3_active = ifelse(cohab3 == 1 & curr3 == 1, 1, 0),
  casual3_active = ifelse(casual3 == 1 & curr3 == 1, 1, 0)
)

# classify one-times as once == 1 and curr == 0 (had sex with once and does not consider to be a current partner)
m <- m |> mutate(
  inst1 = ifelse(once1 == 1 & curr1 == 0, 1, 0),
  inst2 = ifelse(once2 == 1 & curr2 == 0, 1, 0),
  inst3 = ifelse(once3 == 1 & curr3 == 0, 1, 0)
)

# summary current degree
m <- m |>
  rowwise() |>
  mutate(
    deg = sum(curr1, curr2, curr3, na.rm = T),
    deg_mar = sum(mar1_active, mar2_active, mar3_active, na.rm = T),
    deg_cohab = sum(cohab1_active, cohab2_active, cohab3_active, na.rm = T),
    deg_casual = sum(casual1_active, casual2_active, casual3_active, na.rm = T),
    deg_inst = sum(inst1, inst2, inst3, na.rm = T),
    sum_parts = sum(!is.na(curr1), !is.na(curr2), !is.na(curr3))
  )

m <- m |>
  mutate(deg_main = ifelse((deg_mar | deg_cohab) == 1, 1, 0))

m <- m |> mutate(partsdiff = parts1yr - 3) # max # of detailed partners

# assume all partners in last year beyond 3 were one-time
# (applies to very few participants)
m <- m |>
  mutate(deg_inst_high = ifelse(partsdiff > 0, deg_inst + partsdiff, deg_inst)) |>
  mutate(deg_inst_high = ifelse(is.na(deg_inst_high), 0, deg_inst_high))

## Define Age --------- ############################

# Respondent Age
m <- m |> mutate(age = ifelse(ager < 97, ager, NA))

# Partner Race (for current wives and cohabs)
# NO AGE VARS FOR NON-CWP PARTNERS
m <- m |>
  mutate(
    page1 = ifelse(p1relation %in% c(1, 3), cwpage, NA), # current wife or cohab
    page2 = ifelse(p2relation %in% c(1, 3), cwpage, NA),
    page3 = ifelse(p3relation %in% c(1, 3), cwpage, NA)
  )

# Age Diff - males - females
m <- m |> mutate(
  agediff1 = age - page1,
  agediff2 = age - page2,
  agediff3 = age - page3
)

## Define Race / Ethnicity ------ ###################
m <- m |> mutate(race = hisprace2)

# Partner Race (for current partners)
m <- m |>
  mutate(
    prace1 = ifelse(p1relation %in% c(1, 3), cwpnrace, NA), # current wife or cohab
    prace1 = ifelse(is.na(prace1), p1ynrace1, prace1),
    prace2 = ifelse(p2relation %in% c(1, 3), cwpnrace, NA),
    prace2 = ifelse(is.na(prace2), p1ynrace2, prace2),
    prace3 = ifelse(p3relation %in% c(1, 3), cwpnrace, NA),
    prace3 = ifelse(is.na(prace3), p1ynrace3, prace3),
  )

## Condom Use ------- ########

m <- m |> mutate(
  # Overall in last month
  p_cond_month = NA,
  p_cond_month = ifelse(sex4wk >= 1 & !is.na(cond4wk), cond4wk / sex4wk, p_cond_month) # sex >1, cond4wk not missing
  # p_cond_month = ifelse(p_cond_month>=20, NA, p_cond_month), # one person had 5 sex and 100 condom? ignore
  # Condom use at first and last sex with recent partners
)

# Age Group --------- ################
m <- m |>
  mutate(age = ifelse(age == 50, 49.9, age))

m$age_group <- cut(round(m$age), 7)

m <- m |> ungroup()

## Save out ------ ###########

saveRDS(m, here("data", "surveys", "maledat.rds"))

# Combine and save out full dataset -------- #####
d$female <- 1
m$female <- 0
vars_int <- intersect(colnames(d), colnames(m))

d <- d |> select(all_of(vars_int))
m <- m |> select(all_of(vars_int))

all <- rbind(d, m)

all2 <- all |>
  mutate(
    ageF = ifelse(female == 1, age, 0),
    ageM = ifelse(female == 0, age, 0),
    sqrtage = sqrt(age),
    sex1wk = sex4wk / 4
  ) |>
  filter(deg_mar <= 1, deg_cohab <= 1) |> # exclude 1 person with >1 mar, 6 people w/ >1 cohab
  filter(!(deg_mar == 1 & deg_cohab >= 1)) |> # exclude the 4 people who say they have 1 mar and 1 cohab
  mutate(
    ego = row_number(),
    weight = wgt2017_2019,
    weight_in_sample = wgt2017_2019 * (nrow(dat) / sum(wgt2017_2019))
  ) |>
  mutate(condom_use = ifelse(p_cond_month == 0, "Never",
    ifelse(p_cond_month == 1, "Always",
      ifelse(p_cond_month > 0 &
        p_cond_month < 1,
      "Sometimes", NA
      )
    )
  )) |>
  mutate(condom_use2 = ifelse(p_cond_month < 1, "Never",
    ifelse(p_cond_month == 1, "Always", NA)
  )) |>
  mutate(parts1yr_capped = ifelse(parts1yr >= 2, 2, parts1yr)) |>
  mutate(
    race = ifelse(race == 1, "H",
      ifelse(race == 2, "W",
        ifelse(race == 3, "B",
          ifelse(race == 4, "O", NA)
        )
      )
    )
  )

saveRDS(all2, here("data", "nsfg_wide.rds"))

# nolint end
