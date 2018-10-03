#define BOOST_TEST_MODULE "test_aabb"
#include <boost/test/included/unit_test.hpp>
#include <afmize/collision.hpp>
#include <random>

BOOST_AUTO_TEST_CASE(aabb_sphere)
{
    using sphere = afmize::sphere<double>;
    using aabb   = afmize::aabb<double>;
    using point  = mave::vector<double, 3>;

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<double> uni(-100.0, 100.0);

    for(std::size_t i=0; i<1000; ++i)
    {
        const sphere s{uni(mt), point{uni(mt), uni(mt), uni(mt)}};
        const aabb box = afmize::make_aabb(s);

        BOOST_TEST(box.lower[0] == s.center[0] - s.radius,
                   boost::test_tools::tolerance(1e-6));
        BOOST_TEST(box.lower[1] == s.center[1] - s.radius,
                   boost::test_tools::tolerance(1e-6));
        BOOST_TEST(box.lower[2] == s.center[2] - s.radius,
                   boost::test_tools::tolerance(1e-6));

        BOOST_TEST(box.upper[0] == s.center[0] + s.radius,
                   boost::test_tools::tolerance(1e-6));
        BOOST_TEST(box.upper[1] == s.center[1] + s.radius,
                   boost::test_tools::tolerance(1e-6));
        BOOST_TEST(box.upper[2] == s.center[2] + s.radius,
                   boost::test_tools::tolerance(1e-6));
    }
}

BOOST_AUTO_TEST_CASE(merge_aabb)
{
    using sphere = afmize::sphere<double>;
    using aabb   = afmize::aabb<double>;
    using point  = mave::vector<double, 3>;

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<double> uni(-100.0, 100.0);

    for(std::size_t i=0; i<1000; ++i)
    {
        const aabb box1{point{uni(mt), uni(mt), uni(mt)},
                        point{uni(mt), uni(mt), uni(mt)}};
        const aabb box2{point{uni(mt), uni(mt), uni(mt)},
                        point{uni(mt), uni(mt), uni(mt)}};

        const aabb box3 = afmize::merge_aabb(box1, box2);

        BOOST_TEST(box3.lower[0] == std::min(box1.lower[0], box2.lower[0]),
                   boost::test_tools::tolerance(1e-6));
        BOOST_TEST(box3.lower[1] == std::min(box1.lower[1], box2.lower[1]),
                   boost::test_tools::tolerance(1e-6));
        BOOST_TEST(box3.lower[2] == std::min(box1.lower[2], box2.lower[2]),
                   boost::test_tools::tolerance(1e-6));

        BOOST_TEST(box3.upper[0] == std::max(box1.upper[0], box2.upper[0]),
                   boost::test_tools::tolerance(1e-6));
        BOOST_TEST(box3.upper[1] == std::max(box1.upper[1], box2.upper[1]),
                   boost::test_tools::tolerance(1e-6));
        BOOST_TEST(box3.upper[2] == std::max(box1.upper[2], box2.upper[2]),
                   boost::test_tools::tolerance(1e-6));
    }
}
