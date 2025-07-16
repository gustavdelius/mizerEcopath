

test_that("fillDefaultsFromFishBase fills missing values", {
    skip("Skipping FishBase tests as they are currently too slow.")
    df <- data.frame(
        species = c("Hake", "Mackerel"),
        Scientific_name = c("Merluccius merluccius", "Scomber scombrus"),
        w_max = NA,
        stringsAsFactors = FALSE
    )

    filled <- fillDefaultsFromFishBase(df)

    expect_true(all(!is.na(filled$w_max)))
    expect_true("a" %in% names(filled))
    expect_true("age_mat" %in% names(filled))
})


test_that("fillDefaultsFromFishBase respects overwrite argument", {
    skip()
    df <- data.frame(
        species = "Hake",
        Scientific_name = "Merluccius merluccius",
        a = 999,  # fake wrong value
        w_max = NA,
        stringsAsFactors = FALSE
    )

    # Should leave 'a' as 999
    filled_no_overwrite <- fillDefaultsFromFishBase(df, overwrite = FALSE)
    expect_equal(filled_no_overwrite$a, 999)

    # Should replace 'a' with real FishBase value
    filled_with_overwrite <- fillDefaultsFromFishBase(df, overwrite = TRUE)
    expect_true(filled_with_overwrite$a != 999)
})


test_that("fillDefaultsFromFishBase handles species not in FishBase", {
    skip()
    df <- data.frame(
        species = "UnknownFish",
        Scientific_name = "Fictitious imaginarius",
        w_max = NA,
        stringsAsFactors = FALSE
    )

    filled <- fillDefaultsFromFishBase(df, verbose = FALSE)
    expect_true(is.na(filled$w_max))  # Should remain NA
    expect_true("a" %in% names(filled))  # a should still be present
})


test_that("fillDefaultsFromFishBase works with a custom scientific name column", {
    skip()
    df <- data.frame(
        species = "Mackerel",
        sci_name = "Scomber scombrus",
        w_max = NA,
        stringsAsFactors = FALSE
    )

    filled <- fillDefaultsFromFishBase(df, scientific_name_col = "sci_name")
    expect_true(!is.na(filled$w_max))
    expect_true("a" %in% names(filled))
})
