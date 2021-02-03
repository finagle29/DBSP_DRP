import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
from dbsp_drp import manual_aperture

class ManualTracingGUI:
    """
    GUI for manually identifying object traces

    Takes in specs_dict:

    .. code-block::

        specs_dict[target_name] == {
            'spec': spectrum_2d,
            'edges': [all_left_edges, all_right_edges],
            'traces': traces,
            'fwhms': fwhms
        }

    How to use:
    Displays one spectrum at a time
    use left and right arrow keys to go back and forth between spectra
    press m to lay down a (m)anual trace under your cursor
    press d to (d)elete a manual trace near your cursor
    press c, then left-click and drag to set the region to view the (c)ollapsed flux
    """

    def __init__(self, figure, axes, specs_dict):
        self.axes = axes
        self.figure = figure
        self.canvas = figure.canvas

        plt.rcParams['keymap.back'] = ''              # forward / backward keys to enable (Default: left, c, backspace)
        plt.rcParams['keymap.forward'] = ''           # left handed quick navigation (Default: right, v)

        self.canvas.mpl_connect('key_press_event', self.key_press_callback)
        self.canvas.mpl_connect('button_press_event', self.mouse_press_callback)
        self.canvas.mpl_connect('button_release_event', self.mouse_release_callback)
        self.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)

        # init needs to load data in!!!
        self.specs_dict = specs_dict
        self.targets = list(specs_dict)#.keys()
        self.spec_index = 0

        self.manual_dict = {targ: [] for targ in self.targets}
        self.manual_dict = {targ: {'spat_spec': [],
                                    'fwhm': [],
                                    'bgs': [],
                                    'collapse_region': []} for targ in self.targets}

        self._plotted_manual_traces = []

        self.collapse_region_mode = False
        self.drawing_collapse_region = False
        self.collapse_region = [None, None]
        self.collapse_region_artist = None

        # Initial plot!
        self.plot()
        plt.show()

        # TODO: keep collapse region frames
        # TODO: allow adding manual traces in manual_aperture
        self.manual_dict = {k: v for k, v in self.manual_dict.items() if len(v['spat_spec']) > 0}
        #self.manual_dict = {k: {'spat_spec': v, 'fwhm': [], 'bgs': []} for k, v in self.manual_dict.items() if len(v) > 0}

        # Instead of trying to fit the FWHM,
        # let's launch a manual_aperture GUI for each manual trace!
        for target in self.manual_dict:
            targ_dict = self.specs_dict[target]
            spec = targ_dict['spec']
            traces = targ_dict['traces'] if targ_dict['traces'] is not None else []
            fig = plt.figure()
            ax = fig.add_subplot(111)

            man_ap = manual_aperture.ManualApertureGUI(fig, ax, target, spec,
                traces, targ_dict['edges'], targ_dict['fwhms'],
                self.manual_dict[target]['spat_spec'],
                self.manual_dict[target]['collapse_region'])

            self.manual_dict[target]['bgs'] = man_ap.get_background_areas()
            self.manual_dict[target]['spat_spec'] = man_ap.get_manual_traces()
            self.manual_dict[target]['fwhm'] = man_ap.get_fwhms()

        for target in self.manual_dict:
            self.manual_dict[target]['needs_std'] = (self.specs_dict[target]['traces'] is None)

    @property
    def target(self):
        return self.targets[self.spec_index]

    @property
    def spec(self):
        return self.specs_dict[self.target]['spec']

    @property
    def edges(self):
        return self.specs_dict[self.target]['edges']

    @property
    def traces(self):
        return self.specs_dict[self.target]['traces']

    @property
    def mask(self):
        return self.spec.bpmmask == 0

    @property
    def sky_resid(self):
        return (self.spec.sciimg - self.spec.skymodel) * np.sqrt(self.spec.ivarmodel) * self.mask

    def compute_sky_resid(self, spec):
        return (spec.sciimg - spec.skymodel) * np.sqrt(spec.ivarmodel) * (spec.bpmmask == 0)

    def plot(self):
        from dbsp_drp import p200_arm_redux

        self.axes.clear()
        self.canvas.set_window_title(self.target)

        p200_arm_redux.save_one2dspec(self.axes, self.sky_resid, self.edges, self.traces)
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)

        if self.manual_dict[self.target]['collapse_region']:
            self.collapse_region_artist = self.axes.axhspan(self.manual_dict[self.target]['collapse_region'][0],
                self.manual_dict[self.target]['collapse_region'][1], color='purple', alpha=0.2)

        self.canvas.draw()

    def plot_manual_traces(self):
        # remove self._plotted_manual_traces
        for plotted_trace in self._plotted_manual_traces:
#            self.axes.lines.remove(plotted_trace)
            plotted_trace.remove()
        self._plotted_manual_traces = []

        # TODO: use pypeit order to draw manual traces
        for trace in self.manual_dict[self.target]['spat_spec']:
            spat = float(trace.split(':')[0])
            spec = float(trace.split(':')[1])
            line = self.axes.axvline(spat, c='blue')
            self._plotted_manual_traces.append(line)
            point = self.axes.scatter(spat, spec, c='blue', marker='x')
            self._plotted_manual_traces.append(point)
            #line = self.axes.axvline(spat/440, c='blue')
            #self._plotted_manual_traces.append(line)
        self.canvas.draw()

    def key_press_callback(self, event):
        if not event.inaxes:
            return
        key = event.key

        if key == 'right':
            # display next spectrum
            self.spec_index = (self.spec_index + 1) % len(self.targets)
            self.plot()
            self.plot_manual_traces()
        elif key == 'left':
            # display previous spectrum
            self.spec_index = (self.spec_index - 1) % len(self.targets)
            self.plot()
            self.plot_manual_traces()
        elif key == 'm':
            # save mouse cursor position
            self.manual_dict[self.target]['spat_spec'].append(f'{event.xdata}:{event.ydata}')
            self.plot_manual_traces()
        elif key == 'd':
            manual_traces = self.manual_dict[self.target]['spat_spec']

            if manual_traces:
                best_trace_ix = None
                best_spat = None
                for i, trace in enumerate(manual_traces):
                    this_spat = float(trace.split(':')[0])
                    # 5 pixel threshold?
                    if abs(this_spat - event.xdata) < 10:
                        if best_trace_ix is None:
                            best_trace_ix = i
                            best_spat = this_spat
                        else:
                            if (abs(this_spat - event.xdata) < abs(best_spat - event.xdata)):
                                best_trace_ix = i
                                best_spat = this_spat

                if best_trace_ix is not None:
                    manual_traces.pop(best_trace_ix)
                    self.plot_manual_traces()
        elif key == 'c':
            self.collapse_region_mode = True

    def mouse_press_callback(self, event):
        if (self.collapse_region_mode and event.inaxes):
            self.collapse_region[0] = event.ydata
            self.drawing_collapse_region = True

    def motion_notify_callback(self, event):
        if self.drawing_collapse_region:
            self.collapse_region[1] = event.ydata
            self.draw_dragging_collapse_region()

    def mouse_release_callback(self, event):
        if self.drawing_collapse_region:
            self.drawing_collapse_region = False
            self.collapse_region_mode = False
            self.collapse_region[1] = event.ydata
            self.manual_dict[self.target]['collapse_region'] = self.collapse_region.copy()

    def draw_dragging_collapse_region(self):
        if self.collapse_region_artist is not None:
            self.collapse_region_artist.remove()
        self.collapse_region_artist = self.axes.axhspan(self.collapse_region[0],
            self.collapse_region[1], color='purple', alpha=0.2)
        self.canvas.draw()
