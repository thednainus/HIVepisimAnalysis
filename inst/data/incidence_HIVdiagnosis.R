#incidence of HIV diagnosis in MSM in San Diego
# from e-mail from Samantha Tweeten 16 September 2021


years <- c(1980:2021)
frequency <- c(1, NA, 1, 6, 7, 697, 745, 831, 930, 1141, 1132,
              990, 840, 662, 532, 526, 514, 460, 465, 380, 417,
              460, 509, 439, 464, 450, 474,
              465, 451, 480, 439, 394, 375, 381,380, 376, 368,
              310, 286, 234, 145, NA)

incidenceDiag <- data.frame(time = years, frequency = frequency)
