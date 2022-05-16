# Flask

### Revision comparison

La `Revision comparison` ne marche pas parce que le fichier d'entrée et le fichier de sortie n'ont pas les mêmes noms de dimensions:
`180Ma-1x1.nc` a les dimensions `LONGITUDE` et `LATITUDE`
mais les fichiers en sorties par exemple `bathy` ont les dimensions `nav_lon` et `nav_lat`

idem pour le fichier `180MaScoteseKocsis.nc` où ses dim sont `LON1` et `LAT1`

### Heatflow display

Displaying Heatflow does not work in the app because the file is composed of 3 dims (nav_lon, nav_lat and time_counter), it generates an error because it is necessary that the file to display has 2 dims only nav_lon and nav_lat.

---
# Panel

### Lasso select

This tool requires libgeos library (used by *datashader* and *linked_selection*) (install on panel container):
```
apt-get update
apt-get install -y libgeos-dev
```
or add `RUN apt-get update && apt-get install -y libgeos-dev` in panel dockerfile and dockerfile.dev

Lasso select does not work with **dynamic map**

I've got a code like this

```
ls = link_selections.instance(unselected_alpha=0.4)

graphs = hv.DynamicMap(...)
ls(graphs)
```
leading to an error when trying to use the lasso tool:
```
WARNING:param.dynamic_operation: Callable raised "CallbackError("linked_selection aborted because it could not display selection for all elements: ufunc 'invert' not supported for the input types, and the inputs could not be safely coerced to any supported types according to the casting rule ''safe''.")".
```

=> It seems it is a **BUG** ([Check here](https://github.com/holoviz/holoviews/issues/5117))

>**Note**: I downgraded [Panel, Bokeh, Holoviews] to older version (check below section to get explanation) to make it work but it still doesn't work.
>
>However I noticed it works in the Single App Editor so might be intereting to look at

### --autoreload feature

I added a `Multi_Page_WebApp/services/panel/Dockerfile.dev` with an `--autoreload` option to get the automatic refresh of the app when modifying the code.

Originally, when putting this feature in the Dockerfile, this was not working with a `panel: error: unrecognized arguments: --autoreload`. 

This is because Panel used version `0.10.3` released on *18/01/2021* and **autoreload** feature came out on *17/02/2021*.

To fix this I replaced in the `requirements.viz.txt`:
```
panel ~=0.10.0
bokeh < 2.3
holoviews <1.14
```
by 
```
panel ~= 0.12.0
bokeh ~= 2.4.0
holoviews ~= 1.14.0
```

>**Note**: With this newer versions of [Panel, Bokeh, Holoviews] *lasso tool* does not work anymore (check above section). 
>
>So I downgraded [Panel, Bokeh, Holoviews] to older version to make it work but it still doesn't work.
>
>However I noticed it works in the Single App Editor so might be intereting to look at

---

# VS Code for debbuging the app

### Debbuging container named *< container >* using attach to process id

Steps to connect 
- `cmd + shift + p` type `Remote-Containers: attach to Running Container`
- Select a container *< container >*, a new VSCode window will open
- Then go to `Run & Debug` panel 
- Select debbuging configuration in the top left (next to the green arrow)
- Clic on the green arrow to launch debbuging
This leads to a window with error `Timed out waiting for debug server to connect.` 

For the debug config in the `launch.json` I added `"logToFile": true` to get a log file related to debbuging.

Log file is created here `/root/.vscode-server/extensions/ms-python.python-2022.0.1814523869/` in the *< container >* container under files `debugpy.adapter-*.log`

I get this error:
```
I+00000.129: Injector[PID=17] output:
             b"stderr: b'/bin/sh: 1: gdb: not found\\n'"
```
So in the container I run:
```
apt-get update
apt-get install gdb
```
Then I also add to install `lldb` on the mac.

After this I got a new error linked to `ptrace: Operation not permitted`.
So I added in the `docker-compose.dev.yml` to the *< container >* name section:
```
cap_add:
    - SYS_PTRACE
```
... leading to a new error:
```
I+00000.252: Injector[PID=18] output:
             b'stderr: b\'No symbol "DoAttach" in current context.\\n\''
```

**No solution found** ([check here](https://github.com/microsoft/debugpy/issues/762)) 
=> debugpy doesn't support *arm64* architecture

---

# In Calculate Weights step (for mosaix branch): make file selection disapear when MOSAIX selected and enable Start processing button.

`calculate_weights.html`
```
<!--MODIFY-->
<form>
    ...
    <input type="radio" name="engine" id="engineMosaic" value="mosaic" checked> <!--instead of id="engine" for MOSAIC radio box-->
    ...
    <input type="radio" name="engine" id="engineMosaix" value="mosaic"> <!--instead of id="engine" for MOSAIX radio box-->
    ...
    <label for="coordsfile" id="uploadcoordslabel">Coordinates File:</label> <!--instead of <label for="coordsfile">Coordinates File:</label>-->
</form>

...
<script>
    // --------------MODIFY----------------
    $(document).ready(function () {
        function isValid() {
            {% if not has_subbasins %}
                return false
            {% endif %}
            {% if not has_bathy %}
                return false
            {% endif %}
            if (document.getElementById("engineMosaic").checked) {                    // this line is added
                if (document.getElementById("coordsFileUpload").files.length < 1) {
                    return false
            }
            }
            return true
        }
    // ------------------------------------
    ...

    // --------------ADD----------------
    $('#engineMosaic').change(function () {
        updateProcessing();
        document.getElementById("uploadcoords").style.display = "block";
        document.getElementById("uploadcoordslabel").style.display = "block";
    })
    $('#engineMosaix').change(function () {
        updateProcessing();
        document.getElementById("uploadcoords").style.display = "none";
        document.getElementById("uploadcoordslabel").style.display = "none";
    })
    // ---------------------------------
    ...
    </script>
```

---
flask == 2.1.1
simplejson == 3.17.6
blinker == 1.4
python-dotenv
pika == 1.2.0


flask == 2.1.1
simplejson == 3.17.6
blinker == 1.4
python-dotenv == 0.20.0
pika == 1.2.0

flask == 2.1.2
simplejson == 3.17.6
blinker == 1.4
python-dotenv == 0.20.0
pika == 1.2.1