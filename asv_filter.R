library("tidyverse")

asv_filter <- function(asv_read, asv_md, sample_md, min_read_num, min_read_pct, min_sample_num) {
  md <- asv_read |>
    pivot_longer(-SAMPLE_ID, names_to = "ASV_ID", values_to = "READ_NUM") |>
    filter(READ_NUM > 0) |>
    # add ASV metadata
    left_join(asv_md, by = "ASV_ID") |>
    # add sample metadata
    left_join(sample_md, by = "SAMPLE_ID") |>
    # add read percent per sample
    group_by(SAMPLE_ID) |>
    mutate(READ_PCT = READ_NUM / sum(READ_NUM) * 100, .after = READ_NUM) |>
    ungroup() |>
    # add read percent per species per sample
    group_by(SAMPLE_ID) |>
    mutate(READ_PCT_SP = READ_NUM / sum(READ_NUM) * 100, .after = READ_PCT) |>
    ungroup() |>
    # add number of samples an ASV is present in
    group_by(ASV_ID) |>
    mutate(ASV_SAMPLE_N = n(), .after = ASV_ID) |>
    ungroup()
  
  # first filter by min read num or pct
  md_read_filtered <- md |>
    filter(READ_NUM >= min_read_num) |>
    filter(READ_PCT >= min_read_pct)
  
  # calculate number of samples ASV present in after filter
  # and then filter by min sample num
  md_filtered <- md_read_filtered |>
    group_by(ASV_ID) |>
    mutate(ASV_SAMPLE_N_FILTERED = n(), .after = ASV_SAMPLE_N) |>
    ungroup() |>
    filter(ASV_SAMPLE_N_FILTERED >= min_sample_num)
  
  # format filtered asv_read, asv_md for output
  asv_read_filtered <- md_filtered |>
    select(SAMPLE_ID, ASV_ID, READ_NUM) |>
    pivot_wider(names_from = ASV_ID, values_from = READ_NUM, values_fill = 0L)
  
  asv_md_filtered <- asv_md |>
    # add column indicating whether the ASV was filtered
    mutate(ASV_FILTER = ifelse(
      ASV_ID %in% unique(md_filtered$ASV_ID),
      yes = "KEPT", no = "REMOVED"
    ), .after = ASV_ID) |>
    # add number of samples ASV present in before and after filter
    left_join(
      select(md, ASV_ID, ASV_SAMPLE_N) |> distinct(), by = "ASV_ID"
    ) |>
    left_join(
      select(md_filtered, ASV_ID, ASV_SAMPLE_N_FILTERED) |> distinct(), by = "ASV_ID"
    )
  
  list(md = md_filtered, read = asv_read_filtered, asv = asv_md_filtered, sample = sample_md)
}
