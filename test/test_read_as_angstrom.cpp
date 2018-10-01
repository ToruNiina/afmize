#define BOOST_TEST_MODULE "test_read_as_angstrom"
#include <afmize/input_utility.hpp>
// don't include this first to override BOOST_MPL_LIMIT_LIST_SIZE
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE(read_as_angstrom_default)
{
    {
        const toml::value t(10);
        const double v = afmize::read_as_angstrom<double>(t, "test");
        BOOST_TEST(v == 10.0, boost::test_tools::tolerance(1e-6));
    }
    {
        const toml::value t(10.0);
        const double v = afmize::read_as_angstrom<double>(t, "test");
        BOOST_TEST(v == 10.0, boost::test_tools::tolerance(1e-6));
    }
}

BOOST_AUTO_TEST_CASE(read_as_angstrom_pm)
{
    {
        const toml::value t("100.0pm");
        const double v = afmize::read_as_angstrom<double>(t, "test");
        BOOST_TEST(v == 1.0, boost::test_tools::tolerance(1e-6));
    }
    {
        const toml::value t("1_000.000_000pm");
        const double v = afmize::read_as_angstrom<double>(t, "test");
        BOOST_TEST(v == 10.0, boost::test_tools::tolerance(1e-6));
    }
}

BOOST_AUTO_TEST_CASE(read_as_angstrom_angstrom)
{
    {
        const toml::value t("1.0angstrom");
        const double v = afmize::read_as_angstrom<double>(t, "test");
        BOOST_TEST(v == 1.0, boost::test_tools::tolerance(1e-6));
    }
    {
        const toml::value t("1_0.000_000angstrom");
        const double v = afmize::read_as_angstrom<double>(t, "test");
        BOOST_TEST(v == 10.0, boost::test_tools::tolerance(1e-6));
    }

    {
        const toml::value t("1.0Å");
        const double v = afmize::read_as_angstrom<double>(t, "test");
        BOOST_TEST(v == 1.0, boost::test_tools::tolerance(1e-6));
    }
    {
        const toml::value t("1_0.000_000Å");
        const double v = afmize::read_as_angstrom<double>(t, "test");
        BOOST_TEST(v == 10.0, boost::test_tools::tolerance(1e-6));
    }
}

BOOST_AUTO_TEST_CASE(read_as_angstrom_nm)
{
    {
        const toml::value t("1.0nm");
        const double v = afmize::read_as_angstrom<double>(t, "test");
        BOOST_TEST(v == 10.0, boost::test_tools::tolerance(1e-6));
    }
    {
        const toml::value t("1_0.000_000nm");
        const double v = afmize::read_as_angstrom<double>(t, "test");
        BOOST_TEST(v == 100.0, boost::test_tools::tolerance(1e-6));
    }
}

BOOST_AUTO_TEST_CASE(read_as_angstrom_um)
{
    {
        const toml::value t("1e-3um");
        const double v = afmize::read_as_angstrom<double>(t, "test");
        BOOST_TEST(v == 10.0, boost::test_tools::tolerance(1e-6));
    }
    {
        const toml::value t("0.000_1um");
        const double v = afmize::read_as_angstrom<double>(t, "test");
        BOOST_TEST(v == 1.0, boost::test_tools::tolerance(1e-6));
    }

    {
        const toml::value t("1e-3μm");
        const double v = afmize::read_as_angstrom<double>(t, "test");
        BOOST_TEST(v == 10.0, boost::test_tools::tolerance(1e-6));
    }
    {
        const toml::value t("0.000_1μm");
        const double v = afmize::read_as_angstrom<double>(t, "test");
        BOOST_TEST(v == 1.0, boost::test_tools::tolerance(1e-6));
    }
}

BOOST_AUTO_TEST_CASE(read_as_angstrom_mm)
{
    {
        const toml::value t("1e-6mm");
        const double v = afmize::read_as_angstrom<double>(t, "test");
        BOOST_TEST(v == 10.0, boost::test_tools::tolerance(1e-6));
    }
    {
        const toml::value t("0.000_000_1mm");
        const double v = afmize::read_as_angstrom<double>(t, "test");
        BOOST_TEST(v == 1.0, boost::test_tools::tolerance(1e-6));
    }
}
