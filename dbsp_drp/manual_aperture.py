import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
from astropy import stats
from pypeit.spec2dobj import Spec2DObj
from dbsp_drp.gui_helpers import HelpTool

plt.rcParams['toolbar'] = 'toolmanager'

class ManualApertureGUI:
    """
    GUI for manually selecting background and signal for object traces

    How to use:
        - to designate background regions, press (b)ackground and then click and drag
        - left-click and drag blue (manual) traces around
        - right-click and drag to change FWHM of blue (manual) traces
        - press (d)elete to delete background regions under the cursor
    """

    def __init__(self, figure: plt.Figure, axes: plt.Axes, target,
                 spec: Spec2DObj, traces, edges, fitted_fwhms, manual_traces,
                 collapse_region):
        # Takes in figure, axes, spectrum
        self.axes = axes
        self.figure = figure
        self.canvas = figure.canvas

        self.canvas.mpl_connect('key_press_event', self.key_press_callback)
        self.canvas.mpl_connect('button_press_event', self.mouse_press_callback)
        self.canvas.mpl_connect('button_release_event', self.mouse_release_callback)
        self.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        self.canvas.mpl_connect('pick_event', self.pick_callback)

        self.target = target

        # init needs to load data in!!!
        self.spec = spec
        self.traces = traces
        self.fitted_fwhms = fitted_fwhms
        self.edges = edges
        self.nspat = spec.sciimg.shape[0]
        self.manual_traces = manual_traces
        self.fwhms = [4.0] * len(manual_traces)

        self.new_background_area = [None, None]
        self.background_areas = []
        self.background_area_vspans = []
        self.background_selection_mode = False
        self.dragging_bg = False
        self.updating_bg_area = None

        self.dragging_trace = False
        self.manual_trace_lines = []
        self.manual_trace_dragging_index = None
        self.manual_trace_dragging = None


        self.dragging_fwhm = False
        self.fwhm_spans = []
        self.fwhm_dragging_index = None
        self.fwhm_dragging = None

        if collapse_region:
            collapse_region = [int(min(collapse_region)), int(max(collapse_region)) + 1]
        else:
            collapse_region = [0, spec.sciimg.shape[0] - 1]
        self.collapse_mask = np.zeros_like(spec.sciimg, dtype=int)
        self.collapse_mask[collapse_region[0] : collapse_region[1]] = 1
        self.plot_spec = spec.sciimg
        collapsed_spec, _, _ = stats.sigma_clipped_stats(self.plot_spec, mask=(self.spec.bpmmask != 0) | (self.collapse_mask != 1), axis=0, sigma=3.0,
                                        cenfunc='median', stdfunc=lambda x, axis: stats.mad_std(x, axis=axis, ignore_nan=True))
        self.collapsed_spec = collapsed_spec

        self.canvas.manager.toolmanager.remove_tool("help")
        self.canvas.manager.toolmanager.add_tool('Help', HelpTool, helptext=self.__doc__)
        self.canvas.manager.toolbar.add_tool('Help', 'help')

        # Initial plot!
        self.plot()
        plt.show()

        self.normalize_bgs()

    @property
    def mask(self):
        return (self.spec.bpmmask == 0) & (self.collapse_mask == 1)

    def plot(self):
        if len(self.axes.lines) > 0:
            ylim = self.axes.get_ylim()
            xlim = self.axes.get_xlim()
        else:
            xlim = None
            ylim = None

        self.axes.clear()
        self.canvas.set_window_title(self.target)

        masked_spec = np.ma.array(self.plot_spec, mask=self.mask == 0)

        self.axes.plot(np.mean(masked_spec, axis=0), label='mean flux')
        self.axes.plot(np.median(masked_spec, axis=0), label='median flux')


        self.axes.plot(self.collapsed_spec, label='sig clipped flux')

        all_left, all_right = self.edges

        for i in range(all_left.shape[1]):
            self.axes.axvline(np.mean(all_left[:, i]), c='green', lw=1)
            self.axes.axvline(np.mean(all_right[:, i]), c='red', lw=1)

        for i, trace in enumerate(self.traces):
            spat = np.mean(trace)
            self.axes.axvline(spat, c='orange', lw=1)
            self.axes.axvspan(spat - self.fitted_fwhms[i], spat + self.fitted_fwhms[i], color='orange', alpha=0.2)

        self.manual_trace_lines = []
        self.fwhm_spans = []
        for i, trace in enumerate(self.manual_traces):
            spat = float(trace.split(':')[0])
            line = self.axes.axvline(spat, c='blue', lw=1, picker=True)
            line.set_pickradius(5)
            self.manual_trace_lines.append(line)
            self.fwhm_spans.append(self.axes.axvspan(spat - self.fwhms[i], spat + self.fwhms[i], color='blue', alpha=0.2, picker=5))

        self.background_area_vspans = []
        for bg_area in self.background_areas:
            self.background_area_vspans.append(self.axes.axvspan(bg_area[0], bg_area[1], color='gray'))

        bg_patch = mpl.patches.Patch(color='gray', label='background region')
        blue_line = mpl.lines.Line2D([], [], color='blue', lw=1, label='manual trace')
        orange_line = mpl.lines.Line2D([], [], color='orange', lw=1, label="auto-id'd trace")

        handles, _ = self.axes.get_legend_handles_labels()

        self.axes.legend(handles = handles + [bg_patch, blue_line, orange_line])

        self.axes.set_xlabel("Spatial pixel position")
        self.axes.set_ylabel("Counts")


        if xlim is not None:
            self.axes.set_ylim(ylim)
            self.axes.set_xlim(xlim)
        self.canvas.draw()

    def drag_draw_bg(self):
        if self.updating_bg_area is not None:
            self.updating_bg_area.remove()
        if self.new_background_area[0] is not None and self.new_background_area[1] is not None:
            self.updating_bg_area = self.axes.axvspan(self.new_background_area[0], self.new_background_area[1], color='gray')
        self.canvas.draw()

    def drag_draw_manual_trace(self):
        new_xdata = np.ones(2) * float(self.manual_traces[self.manual_trace_dragging_index].split(':')[0])
        self.manual_trace_dragging.set_xdata(new_xdata)

        xy = self.fwhm_spans[self.manual_trace_dragging_index].get_xy()

        xy[(0, 1, 4), 0] = new_xdata[0] - self.fwhms[self.manual_trace_dragging_index]
        xy[(2, 3), 0] = new_xdata[0] + self.fwhms[self.manual_trace_dragging_index]

        xy[(0, 3, 4), 1] = 1
        xy[(1, 2), 1] = 0

        self.fwhm_spans[self.manual_trace_dragging_index].set_xy(xy)

        self.canvas.draw()

    def drag_draw_fwhm(self):
        center = float(self.manual_traces[self.fwhm_dragging_index].split(':')[0])
        xy = self.fwhm_dragging.get_xy()

        xy[(0, 1, 4), 0] = center - self.fwhms[self.fwhm_dragging_index]
        xy[(2, 3), 0] = center + self.fwhms[self.fwhm_dragging_index]

        xy[(0, 3, 4), 1] = 1
        xy[(1, 2), 1] = 0

        self.fwhm_dragging.set_xy(xy)

        self.canvas.draw()

    def key_press_callback(self, event):
        if not event.inaxes:
            return
        key = event.key

        if key == 'b':
            self.background_selection_mode = True
        if key == 'd':
            self.delete(event)

    def mouse_press_callback(self, event):
        if not (event.inaxes and self.background_selection_mode):
            return

        self.new_background_area = [event.xdata, None]
        self.dragging_bg = True

    def mouse_release_callback(self, event):
        if self.background_selection_mode:
            self.new_background_area[1] = event.xdata

            self.background_areas.append(self.new_background_area)

            self.new_background_area = [None, None]

            self.background_selection_mode = False
            self.dragging_bg = False
            self.plot()
        if self.dragging_trace:
            self.dragging_trace = False
            self.manual_trace_dragging = None
            self.plot()
        if self.dragging_fwhm:
            self.dragging_fwhm = False
            self.fwhm_dragging = None
            self.plot()

    def motion_notify_callback(self, event):
        if not (event.inaxes and ((self.background_selection_mode and self.dragging_bg) or self.dragging_trace or self.dragging_fwhm)):
            return

        if self.dragging_bg:
            self.new_background_area[1] = event.xdata
            self.drag_draw_bg()
        if self.dragging_trace:
            self.manual_traces[self.manual_trace_dragging_index] = str(event.xdata) + ':' + self.manual_traces[self.manual_trace_dragging_index].split(':')[1]
            self.drag_draw_manual_trace()
        if self.dragging_fwhm:
            self.fwhms[self.fwhm_dragging_index] = abs(event.xdata - float(self.manual_traces[self.fwhm_dragging_index].split(':')[0]))
            self.drag_draw_fwhm()

    def pick_callback(self, event):
        button = event.mouseevent.button
        if (event.artist in self.manual_trace_lines) and (button == 1):
            self.dragging_trace = True
            self.manual_trace_dragging = event.artist
            self.manual_trace_dragging_index = self.manual_trace_lines.index(event.artist)
        if (event.artist in self.fwhm_spans) and (button == 3):
            self.dragging_fwhm = True
            self.fwhm_dragging = event.artist
            self.fwhm_dragging_index = self.fwhm_spans.index(event.artist)

    def get_background_areas(self):
        bg_arr = []
        all_left, all_right = self.edges
        left_edge = np.median(all_left[0])
        right_edge = np.median(all_right[0])
        slit_width = right_edge - left_edge
        for bg_area in self.background_areas:
            # will need to test
            bg_arr.append(f"{100*(bg_area[0] - left_edge) / slit_width}:{100*(bg_area[1] - left_edge) / slit_width}")
        return ','.join(bg_arr)

    def get_manual_traces(self):
        return self.manual_traces

    def get_fwhms(self):
        return self.fwhms

    def normalize_bgs(self):
        """
        Merges overlapping background areas.

        Call only after the user is done interacting with the GUI
        """
        # first sort each background
        self.background_areas = list(map(sorted, self.background_areas))
        # then sort all backgrounds by first point
        list.sort(self.background_areas, key=lambda bg: bg[0])
        # for each bg area
        i = 1
        while i < len(self.background_areas):
            # if it overlaps with previous area
            if self.background_areas[i][0] < self.background_areas[i-1][1]:
                # edit end of previous area
                self.background_areas[i-1][1] = self.background_areas[i][1]
                # delete this area from list
                self.background_areas.pop(i)
            else:
                # go to next one
                i += 1

    def delete(self, event):
        """
        Delete a background area under (or nearby) the mouse

        The x-tolerance is hardcoded to 5 CCD pixels.
        """
        # loop through background areas
        best_ix = None
        best_diff = 1e99
        for i, bg in enumerate(self.background_areas):
            this_diff = min(abs(bg[0] - event.xdata), abs(bg[1] - event.xdata))
            if (bg[0] <= event.xdata) and (event.xdata <= bg[1]):
                best_diff = 0
                best_ix = i
            elif this_diff < 5:
                if this_diff < best_diff:
                    best_diff = this_diff
                    best_ix = i
        if best_ix is not None:
            self.background_areas.pop(best_ix)
            # remove from plot
            self.background_area_vspans.pop(best_ix).remove()
            self.canvas.draw()
