#define BOOST_TEST_MODULE "test_collision"
#include <boost/test/included/unit_test.hpp>
#include <afmize/collision.hpp>

BOOST_AUTO_TEST_CASE(collision_sphere_sphere)
{
    using sphere = afmize::sphere<double>;
    using point  = mave::vector<double, 3>;

    {
        sphere probe {1.0, point{0.0, 0.0, 10.0}};
        sphere target{1.0, point{0.0, 0.0,  0.0}};

        const auto t = afmize::collision_z(probe, target);
        BOOST_TEST(!std::isnan(t));
        BOOST_TEST(t == -8.0, boost::test_tools::tolerance(1e-6));
    }
    {
        sphere probe {1.0, point{0.0, 0.0, 10.0}};
        sphere target{1.0, point{0.0, 0.0, 20.0}};

        const auto t = afmize::collision_z(probe, target);
        BOOST_TEST(!std::isnan(t));
        BOOST_TEST(t == 8.0, boost::test_tools::tolerance(1e-6));
    }

    {
        sphere probe {1.0, point{0.0, 0.0, 10.0}};
        sphere target{1.0, point{1.0, 1.0,  0.0}};
        const double expect = std::sqrt(2.0) - 10.0;

        const auto t = afmize::collision_z(probe, target);

        BOOST_TEST(!std::isnan(t));
        BOOST_TEST(t == expect, boost::test_tools::tolerance(1e-6));
    }
    {
        sphere probe {1.0, point{0.0, 0.0, 10.0}};
        sphere target{1.0, point{1.0, 1.0, 20.0}};
        const double expect = 10.0 - std::sqrt(2.0);

        const auto t = afmize::collision_z(probe, target);
        BOOST_TEST(!std::isnan(t));
        BOOST_TEST(t == expect, boost::test_tools::tolerance(1e-6));
    }

    {
        sphere probe {1.0, point{ 0.0,  0.0, 10.0}};
        sphere target{1.0, point{10.0, 10.0,  0.0}};

        const auto t = afmize::collision_z(probe, target);
        BOOST_TEST(std::isnan(t));
    }
}
