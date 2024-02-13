import astropy
import astropy.units as u
import calendar
import webbpsf

def get_month_start_end(year, month):
    _, ndays = calendar.monthrange(year, month)
    start_date_str = f"{year:04d}-{month:02d}-01"
    end_date_str = f"{year:04d}-{month:02d}-{ndays:02d}"

    start_date = astropy.time.Time(start_date_str)
    end_date = astropy.time.Time(end_date_str)

    return start_date, end_date


def test_monthly_trending_plot_auto_opdtable():
    trend_table = webbpsf.trending.monthly_trending_plot(2023, 6, instrument='NIRISS', filter='F380M')
    assert(len(trend_table) == 15)


def test_monthly_trending_plot_opdtable_param():
    # Get broad range opdtable to verify our internal filtering works
    start_date, end_date = get_month_start_end(2023, 5)
    start_date2, end_date2 = get_month_start_end(2023, 7)
    # Start a little early, such that we are going to have at least 1 WFS before the start date
    pre_start_date = astropy.time.Time(start_date) - astropy.time.TimeDelta(4 * u.day)

    # Retrieve full OPD table, then trim to the selected time period
    opdtable0 = webbpsf.mast_wss.retrieve_mast_opd_table()
    opdtable0 = webbpsf.mast_wss.deduplicate_opd_table(opdtable0)
    opdtable = webbpsf.mast_wss.filter_opd_table(opdtable0, start_time=pre_start_date, end_time=end_date2)
    trend_table = webbpsf.trending.monthly_trending_plot(2023, 6, opdtable=opdtable, instrument='NIRISS', filter='F380M')
    assert(len(trend_table) == 15)
