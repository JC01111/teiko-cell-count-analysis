-- ============================================================================
-- Loblaw Bio :: cell-count analytics schema (SQLite)
-- ----------------------------------------------------------------------------
-- Design goals
--   1. One row per fact, never per column: cell counts live in a long/EAV-style
--      table so adding a new population (e.g. dendritic cells) is an INSERT into
--      a lookup table, not an ALTER TABLE + rewrite of every query.
--   2. Separate the *who* (subject), the *why* (enrollment: project + condition
--      + treatment + response) and the *what* (sample) so a subject can appear
--      in several studies / treatment arms over time.
--   3. Small dimension tables for controlled vocabularies so new conditions,
--      treatments and sample types are data, not schema.
-- ============================================================================

PRAGMA foreign_keys = ON;

DROP VIEW  IF EXISTS v_analysis_base;
DROP VIEW  IF EXISTS v_sample_frequency;
DROP TABLE IF EXISTS sample_cell_count;
DROP TABLE IF EXISTS sample;
DROP TABLE IF EXISTS enrollment;
DROP TABLE IF EXISTS subject;
DROP TABLE IF EXISTS project;
DROP TABLE IF EXISTS cell_population;
DROP TABLE IF EXISTS sample_type;
DROP TABLE IF EXISTS treatment;
DROP TABLE IF EXISTS condition;

-- ---------------------------------------------------------------- dimensions
CREATE TABLE condition (
    condition_id   INTEGER PRIMARY KEY,
    condition_name TEXT NOT NULL UNIQUE
);

CREATE TABLE treatment (
    treatment_id   INTEGER PRIMARY KEY,
    treatment_name TEXT NOT NULL UNIQUE
);

CREATE TABLE sample_type (
    sample_type_id   INTEGER PRIMARY KEY,
    sample_type_name TEXT NOT NULL UNIQUE   -- PBMC, WB, tumor, ...
);

-- The key to horizontal scalability: populations are rows, not columns.
CREATE TABLE cell_population (
    population_id   INTEGER PRIMARY KEY,
    population_name TEXT NOT NULL UNIQUE,   -- b_cell, cd8_t_cell, ...
    display_name    TEXT,                   -- "B cell"
    sort_order      INTEGER NOT NULL DEFAULT 0
);

-- ------------------------------------------------------------------- studies
CREATE TABLE project (
    project_id   INTEGER PRIMARY KEY,
    project_code TEXT NOT NULL UNIQUE       -- prj1, prj2, ...
);

CREATE TABLE subject (
    subject_id   INTEGER PRIMARY KEY,
    subject_code TEXT NOT NULL UNIQUE,      -- sbj000, ...
    sex          TEXT CHECK (sex IN ('M','F','O','U'))
);

-- One subject enrolled in one project on one treatment arm.
-- Condition / treatment / response / age are properties of THIS enrollment,
-- not of the person forever, so a subject may later join another trial.
CREATE TABLE enrollment (
    enrollment_id INTEGER PRIMARY KEY,
    subject_id    INTEGER NOT NULL REFERENCES subject(subject_id),
    project_id    INTEGER NOT NULL REFERENCES project(project_id),
    condition_id  INTEGER NOT NULL REFERENCES condition(condition_id),
    treatment_id  INTEGER NOT NULL REFERENCES treatment(treatment_id),
    response      TEXT CHECK (response IN ('yes','no') OR response IS NULL),
    age_at_enrollment INTEGER CHECK (age_at_enrollment IS NULL
                                     OR age_at_enrollment BETWEEN 0 AND 130),
    UNIQUE (subject_id, project_id, condition_id, treatment_id)
);

-- ------------------------------------------------------------------- samples
CREATE TABLE sample (
    sample_id      INTEGER PRIMARY KEY,
    sample_code    TEXT NOT NULL UNIQUE,    -- sample00000, ...
    enrollment_id  INTEGER NOT NULL REFERENCES enrollment(enrollment_id),
    sample_type_id INTEGER NOT NULL REFERENCES sample_type(sample_type_id),
    time_from_treatment_start INTEGER       -- days; NULL = untimed / screening
);

-- The fact table. One row per (sample, population).
CREATE TABLE sample_cell_count (
    sample_id     INTEGER NOT NULL REFERENCES sample(sample_id) ON DELETE CASCADE,
    population_id INTEGER NOT NULL REFERENCES cell_population(population_id),
    cell_count    INTEGER NOT NULL CHECK (cell_count >= 0),
    PRIMARY KEY (sample_id, population_id)
) WITHOUT ROWID;

-- --------------------------------------------------------------- index plan
CREATE INDEX idx_enrollment_subject   ON enrollment (subject_id);
CREATE INDEX idx_enrollment_lookup    ON enrollment (condition_id, treatment_id, response);
CREATE INDEX idx_sample_enrollment    ON sample (enrollment_id);
CREATE INDEX idx_sample_type_time     ON sample (sample_type_id, time_from_treatment_start);
CREATE INDEX idx_scc_population       ON sample_cell_count (population_id);

-- ------------------------------------------------------------------- views
-- Part 2 deliverable, expressed once in SQL and reused everywhere.
CREATE VIEW v_sample_frequency AS
SELECT  s.sample_code                       AS sample,
        tot.total_count                     AS total_count,
        cp.population_name                  AS population,
        scc.cell_count                      AS count,
        ROUND(100.0 * scc.cell_count / tot.total_count, 4) AS percentage
FROM        sample_cell_count scc
JOIN        sample            s   ON s.sample_id      = scc.sample_id
JOIN        cell_population   cp  ON cp.population_id = scc.population_id
JOIN       (SELECT sample_id, SUM(cell_count) AS total_count
            FROM   sample_cell_count
            GROUP  BY sample_id) tot ON tot.sample_id = scc.sample_id
WHERE   tot.total_count > 0;

-- Fully denormalised analysis view: one row per sample x population with every
-- clinical attribute attached. This is what the stats + dashboard read.
CREATE VIEW v_analysis_base AS
SELECT  s.sample_code   AS sample,
        sub.subject_code AS subject,
        p.project_code   AS project,
        c.condition_name AS condition,
        t.treatment_name AS treatment,
        e.response       AS response,
        sub.sex          AS sex,
        e.age_at_enrollment AS age,
        st.sample_type_name AS sample_type,
        s.time_from_treatment_start AS time_from_treatment_start,
        cp.population_name AS population,
        scc.cell_count     AS count,
        tot.total_count    AS total_count,
        ROUND(100.0 * scc.cell_count / tot.total_count, 4) AS percentage
FROM        sample_cell_count scc
JOIN        sample          s   ON s.sample_id       = scc.sample_id
JOIN        cell_population cp  ON cp.population_id  = scc.population_id
JOIN        sample_type     st  ON st.sample_type_id = s.sample_type_id
JOIN        enrollment      e   ON e.enrollment_id   = s.enrollment_id
JOIN        subject         sub ON sub.subject_id    = e.subject_id
JOIN        project         p   ON p.project_id      = e.project_id
JOIN        condition       c   ON c.condition_id    = e.condition_id
JOIN        treatment       t   ON t.treatment_id    = e.treatment_id
JOIN       (SELECT sample_id, SUM(cell_count) AS total_count
            FROM   sample_cell_count
            GROUP  BY sample_id) tot ON tot.sample_id = scc.sample_id;
