# Paddock Watch — automatic setup

This builds a live, self-updating map of real paddock conditions for the
property, hosted free on GitHub. Once set up, nobody runs anything: GitHub
rebuilds the data weekly and the family just opens a link.

**What each file is**
- `index.html` — the map the family opens.
- `build_grid.py` — fetches the real DEA satellite data and builds the grid.
- `.github/workflows/build.yml` — tells GitHub to run that weekly.

---

## One-time setup (about 10 minutes, all clicks)

1. **Make a repository.** On github.com, click **New**. Name it e.g.
   `paddock-watch`. Leave it **Public** (free hosting needs public; the link
   is just unlisted, not searchable). Click **Create repository**.

2. **Add the three files.**
   - Click **Add file → Upload files**. Drag in `index.html` and
     `build_grid.py`. Click **Commit changes**.
   - The workflow lives in a folder. Click **Add file → Create new file**.
     In the name box type exactly: `.github/workflows/build.yml`
     Paste the contents of `build.yml` into the editor. Click **Commit changes**.

3. **Turn on Pages (the public link).** Go to **Settings → Pages**. Under
   *Build and deployment*, set **Source** to **Deploy from a branch**, pick
   **Branch: main** and **/ (root)**, and **Save**. After a minute your link
   appears at the top of that page: `https://<you>.github.io/paddock-watch/`.

4. **Build the data the first time.** Go to the **Actions** tab. If it asks,
   click to enable workflows. Pick **Build paddock grid** on the left, then
   **Run workflow → Run workflow**. Wait 3–5 minutes for the green tick.

5. **Open your link.** The grid of real ground cover should now show. It
   refreshes itself every Monday from then on.

---

## Options

- **More lots / bigger property.** Edit `LOTS` in `build.yml`
  (e.g. `"6SP196201,4WLA10"`). Edit `CELL_M` to change cell size in metres
  (1000 = 1 km cells).
- **Rainfall figure.** Go to **Settings → Secrets and variables → Actions →
  New repository secret**. Name it `SILO_EMAIL`, value your email address
  (that's SILO's free API key). Re-run the workflow.
- **Change the schedule.** The `cron` line in `build.yml` is Mondays; change
  it if you like.

## If the first run fails
Open the **Actions** tab, click the failed run, and read the red step. Send
me that error text and I'll fix it. First runs often need one small tweak,
and the log says exactly what.

## Note on privacy
A free GitHub Pages site is public to anyone with the link. For a station
map that's usually fine. A private site needs a paid GitHub plan.
