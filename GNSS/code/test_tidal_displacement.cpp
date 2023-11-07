#include "gtest/gtest.h"

#include "satellite_system/ephemeris/EphemCalculator.hpp"
#include "satellite_system/ephemeris/EphemerisCoefficients.hpp"
#include "satellite_system/reference_system/FrameConverter.hpp"
#include "satellite_system/reference_system/displacement/SolidDisplacementCalculator.hpp"
#include "satellite_system/reference_system/displacement/TidalDisplacementIers2010.hpp"
#include "satellite_system/reference_system/iers2010/ItrsGcrsIERS2010.hpp"
#include "satellite_system/time/TimeConverter.hpp"
#include "satellite_system/utils/Limits.hpp"
#include "satellite_system/utils/MathFunctions.hpp"
#include "third_party/dephem/include/dephem.hpp"

#include "satellite_system/ephemeris/EphemerisCoefficients.hpp"
#include "satellite_system/ephemeris/de405/De405_2010_2040.hpp"
#include "satellite_system/ephemeris/de405/De405_2020_2030.hpp"
#include "satellite_system/ephemeris/de405/De405_2022_2023.hpp"
#include "satellite_system/time/Time.hpp"
#include "satellite_system/time/TimeConverter.hpp"
#include "satellite_system/utils/Utils.hpp"
#include "satellite_system/utils/eop_interpolation/prediction/ZeroDutPredictor.hpp"
#include "satellite_system/utils/eop_interpolation/prediction/ZeroPoleMotion.hpp"
#include "tests/ephemeris/CalcephShell.hpp"
#include "tests/macroses_for_test.hpp"
#include "tests/time_and_frame_test_data/earth_rotation_source.hpp"

namespace ReferenceSystem = SatelliteSystem::ReferenceSystem;
namespace Time = SatelliteSystem::Time;
namespace Math = SatelliteSystem::Utils::Math;
namespace Utils = SatelliteSystem::Utils;

using TimeTdb = Time::Time<SatelliteSystem::Time::TimeScalesEnum::TDB_SCALE>;
using TimeUt1 = Time::Time<SatelliteSystem::Time::TimeScalesEnum::UT1_SCALE>;


using TimeConverter = SatelliteSystem::Time::TimeConverterPrecise;
using DutInterpolator = TimeConverter::DutInterpolator;
using FrameConverter = SatelliteSystem::ReferenceSystem::GcrsToItrsPrecise;
using EopInterpolator = FrameConverter::EopInterpolator;

using TimeTDB = Time::Time<Time::TimeScalesEnum::TDB_SCALE>;
using TimeTT = Time::Time<Time::TimeScalesEnum::TT_SCALE>;
using TimeUT1 = Time::Time<Time::TimeScalesEnum::UT1_SCALE>;

const std::string FILE_PATH = __FILE__;
// Следующая переменная - число символов, которые мы "отрежем" от абсолютного пути текущего файла,
// чтобы получить абсолютный путь к папке satellite_system/
// Для разных файлов это число разное
const int amountOfSymbolsBeforeRootDirectory = 41;
const std::string DIR_PATH = FILE_PATH.substr(0, FILE_PATH.size() - 27 - 36);

using namespace SatelliteSystem;

using TimeScalesEnum = Time::TimeScalesEnum;
using Calendar = Time::Calendar;

constexpr scalar METERS_IN_KILOMETERS = 1000;

class EPHEM_CALCULATOR_TEST : public ::testing::Test {
protected:
    const scalar POSITION_TOLERANCE = 2e4 * Utils::NumericLimits<scalar>::relativeTolerance;
    const scalar VELOCITY_TOLERANCE = 2e4 * Utils::NumericLimits<scalar>::relativeTolerance;

    const indexType dHour = 4;
    const indexType dMinute = 15;

    const Ephemeris::De405_2020_2030 container{};

    std::unique_ptr<Ephemeris::EphemCalculator<Ephemeris::De405_2020_2030>> ephemCalculatorPtr;
    std::unique_ptr<Ephemeris::EphemCalculator<Ephemeris::Calceph>> ephemCalcephPtr;

    using TimeConverter = Time::TimeConverterSimple;
    using DutInterpolator = TimeConverter::DutInterpolator;
    const DutInterpolator data = DutInterpolator::buildDutInterpolator(mjd[0], dut).value();
    const TimeConverter timeConverter = TimeConverter(data);

    void SetUp() override {
        ephemCalculatorPtr = std::make_unique<Ephemeris::EphemCalculator<Ephemeris::De405_2020_2030>>(container);

        ephemCalcephPtr = std::make_unique<Ephemeris::EphemCalculator<Ephemeris::Calceph>>(
            std::move(*Ephemeris::EphemCalculator<Ephemeris::Calceph>::buildEphemCalculator(
                DIR_PATH + "resources/ephemeris/de405.bin")));
    }

    using AllBodiesPack =
        BodiesPack<CelestialBodiesEnum::MERCURY, CelestialBodiesEnum::VENUS, CelestialBodiesEnum::EARTH,
                   CelestialBodiesEnum::MARS, CelestialBodiesEnum::JUPITER, CelestialBodiesEnum::SATURN,
                   CelestialBodiesEnum::URANUS, CelestialBodiesEnum::NEPTUNE, CelestialBodiesEnum::PLUTO,
                   CelestialBodiesEnum::MOON, CelestialBodiesEnum::SUN, CelestialBodiesEnum::SOLAR_SYSTEM_BARYCENTER,
                   CelestialBodiesEnum::EARTH_MOON_BARYCENTER>;

    const std::function<scalar(const Vector3d&, const Vector3d&)> calcError = [](const Vector3d& first,
                                                                                 const Vector3d& second) {
        return first.norm() != 0 ? (first - second).norm() / first.norm() : (first - second).norm();
    };

    template <bool calculateVelocity, typename FirstCalculator, typename SecondCalculator>
    void compareTwoCalculators(const FirstCalculator& first, const SecondCalculator& second, const indexType startYear,
                               const indexType endYear) {
        for (const auto& body : AllCelestialBodiesEnum) {
            for (indexType year = startYear; year < endYear; year++) {
                for (indexType month = 1; month < 13; month++) {
                    for (indexType day = 1; day < 29; day++) {
                        for (indexType hour = 0; hour < 24; hour += dHour) {
                            for (indexType minute = 0; minute < 60; minute += dMinute) {
                                const Calendar calendar = {year, month, day, hour, minute, 0.};
                                const auto timeTdb = *Time::buildFromCalendar<TimeScalesEnum::TDB_SCALE>(calendar);

                                if constexpr (calculateVelocity) {
                                    const auto firstResultExp =
                                        first.calcPositionAndVelocity(timeTdb, body, CelestialBodiesEnum::EARTH);
                                    CREATE_CONST_EXPECTED_OR_ASSERT(firstResultExp, firstResult);

                                    const auto secondResultExp =
                                        second.calcPositionAndVelocity(timeTdb, body, CelestialBodiesEnum::EARTH);
                                    CREATE_CONST_EXPECTED_OR_ASSERT(secondResultExp, secondResult);

                                    const scalar posError = calcError(firstResult.position, secondResult.position);
                                    ASSERT_NEAR(posError, 0, POSITION_TOLERANCE);

                                    const scalar velError = calcError(firstResult.velocity, secondResult.velocity);
                                    ASSERT_NEAR(velError, 0, VELOCITY_TOLERANCE);
                                } else {
                                    const auto firstResultExp =
                                        first.calcPosition(timeTdb, body, CelestialBodiesEnum::EARTH);
                                    CREATE_CONST_EXPECTED_OR_ASSERT(firstResultExp, firstResult);

                                    const auto secondResultExp =
                                        second.calcPosition(timeTdb, body, CelestialBodiesEnum::EARTH);
                                    CREATE_CONST_EXPECTED_OR_ASSERT(secondResultExp, secondResult);

                                    const scalar posError = calcError(firstResult, secondResult);
                                    ASSERT_NEAR(posError, 0, POSITION_TOLERANCE);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    template <bool calculateVelocity, typename Calculator>
    void compareCalculatorWithCalceph(const Calculator& calculator, const indexType startYear,
                                      const indexType endYear) {
        for (const auto& body : AllCelestialBodiesEnum) {
            for (indexType year = startYear; year < endYear; year++) {
                for (indexType month = 1; month < 13; month++) {
                    for (indexType day = 1; day < 29; day++) {
                        for (indexType hour = 0; hour < 24; hour += dHour) {
                            for (indexType minute = 0; minute < 60; minute += dMinute) {
                                const Calendar calendar = {year, month, day, hour, minute, 0.};
                                const auto timeTdb = *Time::buildFromCalendar<TimeScalesEnum::TDB_SCALE>(calendar);

                                const auto timeTt = *timeConverter.convert<Time::TimeScalesEnum::TT_SCALE>(timeTdb);

                                if constexpr (calculateVelocity) {
                                    const auto calcephResultExp = ephemCalcephPtr->calcPositionAndVelocity(
                                        timeTt, timeConverter, body, CelestialBodiesEnum::EARTH);
                                    CREATE_CONST_EXPECTED_OR_ASSERT(calcephResultExp, calcephResult);

                                    const auto containerResultExp =
                                        calculator.calcPositionAndVelocity(timeTdb, body, CelestialBodiesEnum::EARTH);
                                    CREATE_CONST_EXPECTED_OR_ASSERT(containerResultExp, containerResult);

                                    const scalar posError = calcError(containerResult.position, calcephResult.position);
                                    ASSERT_NEAR(posError, 0, POSITION_TOLERANCE);

                                    const scalar velError = calcError(calcephResult.velocity, containerResult.velocity);
                                    ASSERT_NEAR(velError, 0, VELOCITY_TOLERANCE);
                                } else {
                                    const auto calcephResultExp = ephemCalcephPtr->calcPosition(
                                        timeTt, timeConverter, body, CelestialBodiesEnum::EARTH);
                                    CREATE_CONST_EXPECTED_OR_ASSERT(calcephResultExp, calcephResult);

                                    const auto containerResultExp =
                                        calculator.calcPosition(timeTdb, body, CelestialBodiesEnum::EARTH);
                                    CREATE_CONST_EXPECTED_OR_ASSERT(containerResultExp, containerResult);

                                    const scalar posError = calcError(containerResult, calcephResult);
                                    ASSERT_NEAR(posError, 0, POSITION_TOLERANCE);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
};

TEST_F(EPHEM_CALCULATOR_TEST, MANY_TIMES) {
    // 2022-11-20 08:00:30.0
    // 2022-11-20 23:59:30.0
    const std::string FILE_PATH = __FILE__;
    std::string DIR_PATH = FILE_PATH.substr(0, FILE_PATH.size() - 27 - 36);
    // /Users/arseniy/work/satellite_system/satellite_system/resources/ephemeris/de405.bin
    auto desktopContainer = Ephemeris::EphemerisCoefficients::build(DIR_PATH + "resources/ephemeris/de405.bin");
    ASSERT_TRUE(desktopContainer);

    const auto desktopEphemCalculator =
        Ephemeris::EphemCalculator<Ephemeris::EphemerisCoefficients>(std::move(*desktopContainer));

    // 1 Frame and time data
    using Dut = SatelliteSystem::Utils::EopInterpolation::ZeroDutPredictor;
    using Dtr = SatelliteSystem::Time::ApproxDtrCalculator;
    using DutCorr = SatelliteSystem::Time::DutNoCorrections;
    using TimeConverter = SatelliteSystem::Time::TimeConverterIers2010<Dut, DutCorr, Dtr>;
    TimeConverter timeConverter = TimeConverter(Dut{});

    using PoleMotionCorr = SatelliteSystem::ReferenceSystem::PoleMotionNoCorrections;
    using PoleMotion = SatelliteSystem::Utils::EopInterpolation::ZeroPoleMotionPredictor;
    using FrameConverter = SatelliteSystem::ReferenceSystem::ItrsGcrsIERS2010<PoleMotion, PoleMotionCorr>;
    FrameConverter frameConverter = FrameConverter(PoleMotion{});

    SatelliteSystem::Vector3d displacement;

    using DispCalculator = ReferenceSystem::TidalDisplacementCalculatorIers<ReferenceSystem::TidalDisplacementIers2010>;
//        const SatelliteSystem::Vector3d pointItrs(0, 0, 6356e3); // polar
    SatelliteSystem::Vector3d pointItrs(6378e3, 0, 0);  // equator
    const SatelliteSystem::Vector3d moonPositionItrs(-312468450.131567, 0, -312468450.131567);
    const SatelliteSystem::Vector3d sunPositionItrs(0, 0, 137859926952.015);

    const auto timeUt1Start =
        Time::buildFromCalendar<SatelliteSystem::Time::TimeScalesEnum::UT1_SCALE>(2022, 11, 1, 1, 0, 0).value();

    std::fstream file;
    DIR_PATH = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    file.open(DIR_PATH + "GNSS.txt", std::ios::out);
    for(int d = 1; d < 25; ++d) {
        for (double h = 1; h < 24; ++h) {
            for (double m = 0; m < 60; ++m) {
                for (double s = 0; s < 60; s += 30) {
                    // конвертируешь время tdb
                    const auto timeUt1Exp = Time::buildFromCalendar<SatelliteSystem::Time::TimeScalesEnum::UT1_SCALE>(
                        2022, 11, d, h, m, s);
                    ASSERT_TRUE(timeUt1Exp);
                    const auto timeUt1 = timeUt1Exp.value();
                    const auto timeTdb = TimeTdb(timeUt1.jdDay(), timeUt1.jdDayPart()).addSeconds(32.183999999999997);
                    const auto timeTT = timeConverter.convert<SatelliteSystem::Time::TimeScalesEnum::TT_SCALE>(timeUt1);
                    // находимшь кватернион между системами
                    const auto moonPositionGcrs = desktopEphemCalculator
                                                      .calcPosition(timeTdb, SatelliteSystem::CelestialBodiesEnum::MOON,
                                                                    SatelliteSystem::CelestialBodiesEnum::EARTH)
                                                      .value();
                    const auto sunPositionGcrs = desktopEphemCalculator
                                                     .calcPosition(timeTdb, SatelliteSystem::CelestialBodiesEnum::SUN,
                                                                   SatelliteSystem::CelestialBodiesEnum::EARTH)
                                                     .value();

                    // конвертируешь вектора с помощью кватерниона
                    Quaterniond quat = frameConverter.eciToEcef(timeTT.value(), timeConverter).value();

                    displacement = DispCalculator::calcDisplacementTideFree(
                        pointItrs, quat._transformVector(moonPositionGcrs), quat._transformVector(sunPositionItrs),
                        timeTdb, timeUt1);

                    file << displacement(0) << " " << displacement(1) << " " << displacement(2) << " "
                         << displacement.norm() << " " << timeUt1 - timeUt1Start << "\n";
                }
            }
        }
    }
    file.close();
}
