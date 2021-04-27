from dbsp_drp.instruments.instrument import Instrument
from dbsp_drp.instruments.dbsp import DBSP
from dbsp_drp.instruments.utils import guess_instrument_from_file

instruments = Instrument.__subclasses__()
