-- include-code.lua: fill code blocks on the slides from files in the repo.
--
-- A code block carrying `include="path"` has its body replaced by that file.
-- With `snippet="name"` as well, only the lines between a `start snippet name`
-- and an `end snippet name` marker line are used (the markers themselves are
-- dropped), so a runnable script can carry plumbing that stays off the slide:
--
--     ```{.r include="scripts/03_least_squares.R" snippet="ls_fitting"}
--     ```
--
-- `start-line="n"` and `end-line="m"` (1-based, inclusive) keep only that
-- range of the file or snippet, for splitting a long listing over slides.
--
-- Same attributes as quarto-ext/include-code-files, with three differences
-- that matter for a deck: a missing file or snippet stops the render instead
-- of silently dumping the whole file onto a slide, marker names are matched
-- exactly (so `ls_fitting` never picks up `ls_fitting_2`), and line numbers on
-- the slide start at 1 so `code-line-numbers` highlights line up with them.

local function read_lines(path)
  local fh = io.open(path)
  if not fh then
    error("include-code: cannot open " .. path)
  end
  local lines = {}
  for line in fh:lines() do
    lines[#lines + 1] = line
  end
  fh:close()
  return lines
end

local function snippet(lines, name, path)
  local out, inside = {}, false
  for _, line in ipairs(lines) do
    if inside then
      if line:match("end snippet%s+(%S+)%s*$") == name then
        return out
      end
      out[#out + 1] = line
    elseif line:match("start snippet%s+(%S+)%s*$") == name then
      inside = true
    end
  end
  error(("include-code: no complete snippet '%s' in %s"):format(name, path))
end

function CodeBlock(cb)
  local path = cb.attributes.include
  if not path then
    return nil
  end
  local lines = read_lines(path)
  if cb.attributes.snippet then
    lines = snippet(lines, cb.attributes.snippet, path)
  end
  if cb.attributes["start-line"] or cb.attributes["end-line"] then
    local first = tonumber(cb.attributes["start-line"]) or 1
    local last = tonumber(cb.attributes["end-line"]) or #lines
    local kept = {}
    for i = first, math.min(last, #lines) do kept[#kept + 1] = lines[i] end
    lines = kept
  end
  cb.attributes.include = nil
  cb.attributes.snippet = nil
  cb.attributes["start-line"] = nil
  cb.attributes["end-line"] = nil
  cb.text = table.concat(lines, "\n")
  return cb
end
