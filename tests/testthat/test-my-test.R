test_that("multiplication works", {

    func1 <- function() {
        return( invisible(4) )
    }

    expect_equal( func1(), 4 )



    func2 <- function() {

        print("abc")
        print("abadasdc")

    }

    expect_output( func2(), "abc" )
})
