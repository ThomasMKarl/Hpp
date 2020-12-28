#include "gtest/gtest.h"
#include "grid/grid.h"

namespace HB {
// The fixture for testing class Foo.
class GridTest : public ::testing::Test
{
 protected:
  // You can remove any or all of the following functions if their bodies would
  // be empty.

  GridTest()
  {
     // You can do set-up work for each test here.
  }

  ~GridTest() override
  {
     // You can do clean-up work that doesn't throw exceptions here.
  }

  // If the constructor and destructor are not enough for setting up
  // and cleaning up each test, you can define the following methods:

  void SetUp() override
  {
     // Code here will be called immediately after the constructor (right
     // before each test).
  }

  void TearDown() override
  {
     // Code here will be called immediately after each test (right
     // before the destructor).
  }

  // Class members declared here can be used by all tests in the test suite
  // for Foo.
};

TEST_F(GridTest, calcNeighbourTable_NeighbourIndicesAreCorrect)
{
  EXPECT_EQ(0, 0);
}

TEST_F(GridTest, hotStart_AngleIsValid)
{
  EXPECT_EQ(0, 0);
}

TEST_F(GridTest, hotStart_IsingDirectionIsValid)
{
  EXPECT_EQ(0, 0);
}

TEST_F(GridTest, hotStart_SpinsAreValid)
{
  EXPECT_EQ(0, 0);
}

TEST_F(GridTest, coldStart_AngleIsDown)
{
  EXPECT_EQ(0, 0);
}

TEST_F(GridTest, coldStart_IsingDirectionIsDown)
{
  EXPECT_EQ(0, 0);
}

TEST_F(GridTest, coldStart_SpinIsDown)
{
  EXPECT_EQ(0, 0);
}

}
